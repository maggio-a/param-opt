#ifndef PARALLEL_H
#define PARALLEL_H

#include <vector>
#include <memory>

#include <mutex>
#include <deque>
#include <condition_variable>
#include <atomic>
#include <thread>

#include "mesh.h"
#include "mesh_graph.h"
#include "parameterization.h"

#include <vcg/complex/algorithms/hole.h>

namespace parallel {

inline void Init()
{
    vcg::tri::MinimumWeightEar<Mesh>::ThreadSafeInit();
}

template <typename T>
class ConcurrentQueue {

private:

    const std::size_t _capacity; // 0 means unbounded
    std::mutex _mtx;
    std::condition_variable _resume;
    std::deque<T> _deque;

public:

    ConcurrentQueue(std::size_t capacity=0) : _mtx{}, _deque{}, _capacity{capacity} { }
    ConcurrentQueue(const ConcurrentQueue<T>& other) = delete;
    ConcurrentQueue<T>& operator=(const ConcurrentQueue<T>& other) = delete;

    bool IsEmpty()
    {
        std::lock_guard<std::mutex> lock{_mtx};
        return _deque.empty();
    }

    void Put(const T& v)
    {
        std::unique_lock<std::mutex> lock{_mtx};
        while (_capacity > 0 && _deque.size() > _capacity)
            _resume.wait(lock);
        _deque.push_back(v);
    }

    bool Get(T& v)
    {

        std::lock_guard<std::mutex> lock{_mtx};
        if (!_deque.empty()) {
            v = _deque.front();
            _deque.pop_front();
            if (_capacity > 0 && _deque.size() < 0.25*_capacity)
                _resume.notify_all();
            return true;
        } else return false;
    }

};

class WorkerPool {

    struct TaskType {
        ChartHandle chart;
        bool warmStart;
    };

    using TaskQueue = ConcurrentQueue<TaskType>;

    int threadnum;

    GraphManager& graphMgr; // should be a unique_ptr
    std::mutex gmMutex;

    std::unique_ptr<TaskQueue> globalQueue;
    std::vector<std::unique_ptr<TaskQueue>> workerQueue;

    std::atomic<int> numActiveTasks; // active tasks are either tasks in a queue or tasks being currently processed

    std::vector<std::thread> threads;

    std::atomic<int> numFinalCharts;
    std::atomic<int> numBacktrackFailed;
    std::atomic<int> numBacktrackNonInjective;

    void ParameterizationWorker(int id, ParameterizationStrategy strategy, double injectivityTolerance);
    bool Steal(TaskType& tv, int id);

public:

    WorkerPool(int tn, GraphManager& gm);

    void Run(ParameterizationStrategy strategy, double injectivityTolerance);

};



} // namespace

#endif // PARALLEL_H

