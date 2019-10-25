import multiprocessing as mp

'''
A simple parallel function runner for jupyter notebook
'''

def runas(newname):
    def decorator(f):
        f.__name__ = newname
        return f
    return decorator

def parallel_run(f, data, args=None, kwargs=None):
    '''
    Runs function f over data provided in data in parallel.

    f: function to be run that needs to accept queue as a keyword argument
    data: an iterable producing positional arguments to f, dicts are interpreted
    as kwargs for f.
    args: possible additional positional arguments to be passed to f
    kwargs: keyword arguments to be passed to f
    '''
    output = mp.Queue()
    pool = []
    if kwargs is None:
        kwargs = {'queue': output,}
    else:
        kwargs['queue'] = output
    for i,d in enumerate(data):
        if isinstance(d, dict):
            kwargs.update(d)
        else:
            _args = list(d)
            if args is not None:
                _args.extend(args)
        
        pool.append(mp.Process(target=f, args=_args, kwargs=kwargs))
    
    print("started {} processes ...".format(i+1))
    
    for p in pool:
        p.start()

    for p in pool:
        p.join()

    stat = [output.get() for p in pool]
    return stat
