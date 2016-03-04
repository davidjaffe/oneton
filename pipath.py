#!/usr/bin/env python
'''
generate platform-independent path given posix style path
20160304
'''
import os

class pipath():
    def __init__(self):
        return
    def fix(self,posixpath):
        '''
        return platform-indepent path give posix style path
        handle single string or list of strings as input
        from 
        http://stackoverflow.com/questions/13162372/using-absolute-unix-paths-in-windows-with-python?rq=1
        '''
        if isinstance(posixpath,list):
            fixed = [os.path.abspath(x) for x in posixpath]
        else:
            fixed = os.path.abspath(posixpath)
        return fixed
    def oldfix(self,posixpath):
        '''
        NOT SO GOOD CUZ INITIAL '/' IN posixpath IS STRIPPED OFF
        return platform-indepent path give posix style path
        handle single string or list of strings as input
        '''
        if isinstance(posixpath,list):
            fixed = []
            for x in posixpath:
                s = x.split('/')
                fixed.append( os.path.join(*s))
        else:
            s = posixpath.split('/')
            fixed = os.path.join(*s)
        return fixed
        
