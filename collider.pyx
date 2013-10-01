'''
Collider
======

The collider module contains classes which can be used to test membership
of a point in some space. See individual class documentation for details.

To use it you first have to cython the code. To do that, cd to the directory
containing collider.pyx and::

    python setup.py build_ext --inplace
'''

cimport cython
from libc.stdlib cimport malloc, free
cdef extern from "math.h":
    double round(double val)
    double floor(double val)
    double ceil(double val)

__all__ = ('Collide2DPoly', )


cdef class Collide2DPoly(object):
    ''' Collide2DPoly checks whether a point is within a polygon defined by a
    list of corner points.

    Based on http://alienryderflex.com/polygon/

    For example, a simple triangle::
    >>> collider = Collide2DPoly([10., 10., 20., 30., 30., 10.],\
                                 pre_compute=True)
    >>> (0.0, 0.0) in collider
    False
    >>> (20.0, 20.0) in collider
    True

    The constructor takes a list of x,y points in the form of [x1,y1,x2,y2...]
    as the points argument. These points define the corners of the
    polygon. The boundary is linearly interpolated between each set of points.
    The x, and y values must be floating points.
    The cache argument, if True, will calculate membership for all the points
    so when collide_point is called it'll just be a table lookup.
    '''

    cdef double *cpoints
    cdef double *cconstant
    cdef double *cmultiple
    cdef char *cspace
    cdef double min_x
    cdef double max_x
    cdef double min_y
    cdef double max_y
    cdef int width
    cdef int count

    @cython.cdivision(True)
    def __cinit__(self, points, cache=False, **kwargs):
        cdef int length = len(points)
        if length % 2:
            raise IndexError()
        if length < 4:
            self.cpoints = NULL
            return

        cdef int count = length / 2
        self.count = count
        self.cpoints = <double *>malloc(length * cython.sizeof(double))
        self.cconstant = <double *>malloc(count * cython.sizeof(double))
        self.cmultiple = <double *>malloc(count * cython.sizeof(double))
        cdef double *cpoints = self.cpoints
        cdef double *cconstant = self.cconstant
        cdef double *cmultiple = self.cmultiple
        self.cspace = NULL
        if cpoints is NULL or cconstant is NULL or cmultiple is NULL:
            raise MemoryError()

        self.min_x = min(points[0::2])
        self.max_x = max(points[0::2])
        self.min_y = min(points[1::2])
        self.max_y = max(points[1::2])
        cdef double min_x = floor(self.min_x), min_y = floor(self.min_y)
        cdef int i_x, i_y, j_x, j_y, i
        cdef int j = count - 1, odd, width, height, x, y
        for i in range(length):
            cpoints[i] = points[i]
        if cache:
            for i in range(count):
                cpoints[2 * i] -= min_x
                cpoints[2 * i + 1] -= min_y

        for i in range(count):
            i_x = i * 2
            i_y = i_x + 1
            j_x = j * 2
            j_y = j_x + 1
            if cpoints[j_y] == cpoints[i_y]:
                cconstant[i] = cpoints[i_x]
                cmultiple[i] = 0.
            else:
                cconstant[i] = (cpoints[i_x] - cpoints[i_y] * cpoints[j_x] /
                               (cpoints[j_y] - cpoints[i_y]) +
                               cpoints[i_y] * cpoints[i_x] /
                               (cpoints[j_y] - cpoints[i_y]))
                cmultiple[i] = ((cpoints[j_x] - cpoints[i_x]) /
                               (cpoints[j_y] - cpoints[i_y]))
            j = i
        if cache:
            width = int(ceil(self.max_x) - min_x + 1.)
            self.width = width
            height = int(ceil(self.max_y) - min_y + 1.)
            self.cspace = <char *>malloc(height * width * cython.sizeof(char))
            if self.cspace is NULL:
                raise MemoryError()
            for y in range(height):
                for x in range(width):
                    j = count - 1
                    odd = 0
                    for i from 0 <= i < count:
                        i_y = i * 2 + 1
                        j_y = j * 2 + 1
                        if (cpoints[i_y] < y and cpoints[j_y] >= y or
                            cpoints[j_y] < y and cpoints[i_y] >= y):
                            odd ^= y * cmultiple[i] + cconstant[i] < x
                        j = i
                    self.cspace[y * width + x] = odd
    
    def __dealloc__(self):
        free(self.cpoints)
        free(self.cconstant)
        free(self.cmultiple)
        free(self.cspace)

    @cython.cdivision(True)
    def collide_point(Collide2DPoly self, float x, float y):
        points = self.cpoints
        if points is NULL or not (self.min_x <= x <= self.max_x and
                                  self.min_y <= y <= self.max_y):
            return False
        if self.cspace is not NULL:
            y -= floor(self.min_y)
            x -= floor(self.min_x)
            return self.cspace[int(y) * self.width + int(x)]

        cdef int j = self.count - 1, odd = 0, i, i_y, j_y, i_x, j_x
        for i in range(self.count):
            i_y = i * 2 + 1
            j_y = j * 2 + 1
            if (points[i_y] < y and points[j_y] >= y or
                points[j_y] < y and points[i_y] >= y):
                odd ^= y * self.cmultiple[i] + self.cconstant[i] < x
            j = i
        return odd

    def __contains__(self, point):
        return self.collide_point(*point)
