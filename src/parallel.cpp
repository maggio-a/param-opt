#include "parallel.h"
#include "utils.h"
#include "timer.h"
#include "texture_optimization.h"
#include "texture_rendering.h"

namespace parallel {

WorkerPool::WorkerPool(int tn, GraphManager& gm)
    : threadnum{tn},
      graphMgr{gm},
      gmMutex{},
      numActiveTasks{0},
      numFinalCharts{0},
      numBacktrackFailed{0},
      numBacktrackNonInjective{0}
{
    globalQueue = std::unique_ptr<TaskQueue>(new TaskQueue{});
    for (int i = 0; i < threadnum; ++i) {
        workerQueue.emplace_back(std::unique_ptr<TaskQueue>(new TaskQueue{}));
    }
}

void WorkerPool::Run(ParameterizationStrategy strategy, double injectivityTolerance)
{
    // fill the global queue, sorting tasks by decreasing chart size
    {
        std::vector<TaskType> tasks;
        for (auto entry : graphMgr.Graph()->charts) {
            ChartHandle chart = entry.second;
            Mesh probe;
            MeshFromFacePointers(chart->fpVec, probe);
            if (Parameterizable(probe)) {
                tasks.push_back(TaskType{chart, false});
                numActiveTasks++;
            } else {
                numFinalCharts++;
            }
        }

        std::sort(tasks.begin(), tasks.end(), [](const TaskType& t1, const TaskType& t2) { return t1.chart->FN() > t2.chart->FN(); });
        for (auto& t : tasks)
            globalQueue->Put(t);
    }

    // start the threads
    for (int i = 0; i < threadnum; ++i) {
        threads.emplace_back(std::thread{&WorkerPool::ParameterizationWorker, this, i, strategy, injectivityTolerance});
    }

    // join the threads
    for (std::size_t i = 0; i < threads.size(); ++i) {
        if (threads[i].joinable())
            threads[i].join();
    }

}

bool WorkerPool::Steal(TaskType& tv, int id)
{
    for (int i = 0; i < threadnum; ++i) {
        int index = (id+i+1) % threadnum;
        if (workerQueue[index]->Get(tv))
            return true;
    }
    return false;
}

void WorkerPool::ParameterizationWorker(int id, ParameterizationStrategy strategy, double injectivityTolerance)
{
    bool injectivityCheckRequired;
    if (strategy.scaffold || injectivityTolerance < 0)
        injectivityCheckRequired = false;
    else {
        ensure_condition(injectivityTolerance >= 0 && injectivityTolerance <= 1);
        injectivityCheckRequired = true;
    }

    Timer timer;

    ensure_condition(HasWedgeTexCoordStorageAttribute(graphMgr.Graph()->mesh));
    auto wtcsattr = GetWedgeTexCoordStorageAttribute(graphMgr.Graph()->mesh);

    TaskType task;

    int noTaskIter = 0;
    while (true) {
        if (workerQueue[id]->Get(task) || globalQueue->Get(task) || Steal(task, id)) {

            noTaskIter = 0;

            ChartHandle chart = task.chart;

            //std::cout << "Chart " << chart->id << " - FN=" << chart->FN() << ", FI=" << tri::Index(graph->mesh, chart->Fp()) << std::endl;

            ParameterizationStrategy strat = strategy;
            if (task.warmStart)
                strat.warmStart = true;

            bool parameterized = false;

            Timer t;

            ParameterizerObject po{chart, strat};
            po.Initialize();
            if (po.GetStatus() == ParameterizerObject::Status::Initialized)
                parameterized = po.Parameterize();

            if (parameterized) {
                po.SyncChart();

                bool backtrack = false;
                if (injectivityCheckRequired) {
                    RasterizedParameterizationStats stats = GetRasterizationStats(chart, 1024, 1024);
                    double fraction = stats.lostFragments / (double) stats.totalFragments;
                    if (fraction > injectivityTolerance) {
                        //std::cout << "WARNING: REGION " << chart->id << " HAS OVERLAPS IN THE PARAMETERIZATION (overlap fraction = " << fraction << ")" << std::endl;
                        backtrack = true;
                    }

                }
                if (backtrack) {
                    // split the chart in the mesh graph using the graph manager, and readd
                    // the new task with the warmStart flag set to true

                    std::lock_guard<std::mutex> gmLock{gmMutex};

                    std::vector<ChartHandle> splitCharts;
                    graphMgr.Split(chart->id, splitCharts);
                    std::vector<ChartHandle> newCharts;
                    RecoverFromSplit(splitCharts, graphMgr, newCharts, true);
                    for (auto& c : newCharts) {
                        numActiveTasks++;
                        workerQueue[id]->Put(TaskType{c, strategy.warmStart});
                    }
                    numBacktrackNonInjective++;
                } else {
                    double oldUvArea = chart->OriginalAreaUV();
                    double newUvArea = chart->AreaUV();
                    double scale = std::sqrt(oldUvArea / newUvArea);
                    ensure_condition(scale > 0);
                    vcg::Box2d uvBox = chart->UVBox();
                    for (auto fptr : chart->fpVec) {
                        for (int i = 0; i < 3; ++i) {
                            fptr->WT(i).P() = (fptr->WT(i).P() - uvBox.min) * scale;
                        }
                    }
                    chart->numMerges = 0;
                    numFinalCharts++;
                }
            } else {
                // split the aggregate and restore the original uvs

                std::lock_guard<std::mutex> gmLock{gmMutex};

                bool recover = (chart->numMerges > 0);
                std::vector<ChartHandle> splitCharts;
                graphMgr.Split(chart->id, splitCharts);
                if (recover) {
                    std::vector<ChartHandle> newCharts;
                    RecoverFromFailedInit(splitCharts, graphMgr, newCharts);
                    for (auto& c : newCharts) {
                        numActiveTasks++;
                        workerQueue[id]->Put(TaskType{c, false});
                    }
                } else {
                    //std::cout << "Chart cannot be split, restoring uvs..." << std::endl;
                    for (auto split : splitCharts) {
                        for (auto fptr : split->fpVec) {
                            TexCoordStorage tcs = wtcsattr[fptr];
                            for (int i = 0; i < 3; ++i) {
                                fptr->WT(i) = tcs.tc[i];
                                fptr->V(i)->T() = tcs.tc[i];
                            }
                        }
                    }
                }
            }

            numActiveTasks--;

            std::cout << "Thread " << id << ": task took " << t.TimeElapsed() << " seconds" << std::endl;

        } else {
            noTaskIter++;
            if (numActiveTasks > 0)
                //std::this_thread::yield();
                std::this_thread::sleep_for(std::chrono::milliseconds(10) * std::min(noTaskIter, 10));
            else
                return;
        }
    }
}

} // namespace
