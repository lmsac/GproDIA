import os
import re
import itertools
import pickle
import json

def list_files(path='.', pattern=None, recursive=False, include_dirs=False):
    if recursive:
        return list(itertools.chain.from_iterable(
            [
                os.path.join(t[0], x) 
                for x in (t[1] + t[2] if include_dirs else t[2])
                if pattern is None or \
                re.search(pattern, x) is not None
            ]
            for t in os.walk(path)
        ))
    else:
        return [
             t
             for t in os.listdir(path)
             if pattern is None or \
             re.search(pattern, t) is not None
        ]


def save_json(data, file, **kwargs):
    with open(file, 'w') as f:
        json.dump(data, f, **kwargs)        

def load_json(file, **kwargs):
    with open(file, 'r') as f:
        return json.load(f, **kwargs)
    

def save_pickle(data, file, **kwargs):
    with open(file, 'wb') as f:
        pickle.dump(data, f, **kwargs)
        
def load_pickle(file, **kwargs):
    with open(file, 'rb') as f:
        return pickle.load(f, **kwargs)
    
    