import multiprocessing as mp

'''
A simple parallel function runner for jupyter notebook
'''

def runas(newname):
    def decorator(f):
        f.__name__ = newname
        return f
    return decorator

@runas('main')
def parallel_run(f, data, args=None, kwargs=None):
    '''
    Runs function f over data provided in data in parallel.

    f: function to be run
    data: an iterable producing positional arguments to f
    args: possible additional positional arguments to be passed to f
    kwargs: keyword arguments to be passed to f
    '''
    if __name__ == '__main__':
        output = mp.Queue()
        pool = []
        for i,d in enumerate(data):
            _args = list(d)
            if args is not None:
                _args.extend(args)
            
            if kwargs is not None:
                kwargs = {'queue': output,}
            else:
                kwargs['queue'] = output
            pool.append(mp.Process(target=f, args=_args, kwargs=kwargs))
        
        print("started {} processes ...".format(i+1))
        
        for p in pool:
            p.start()

        for p in pool:
            p.join()

        stat = [output.get() for p in pool]
        return stat
