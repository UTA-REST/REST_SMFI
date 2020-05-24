

def Hist_maker(double r_cut, int Yoff, int Xoff, int shape, ReducedImage):
    cdef int r
    cdef int c
    HIST_VAL = []
    for r in range(shape):
        for c in range(shape):
            rr = (r-Yoff)**2+(c-Xoff)**2
            if rr<=r_cut:
                HIST_VAL.append(ReducedImage[r,c])
    return HIST_VAL