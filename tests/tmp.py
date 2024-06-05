def wcov(v1, v2, w):
    r = np.sum(w * (v1 - np.average(v1, weights=w)) * (v2 - np.average(v2, weights=w))) / np.sum(w)
    return r
#return np.average((v1 - np.average(v1, weights=w)) * (v2 - np.average(v2, weights=w)), weights=w)

def wcorr(v1, v2, w, get_markers=False):
    res = wcov(v1, v2, w) / np.sqrt(wcov(v1, v1, w) * wcov(v2, v2, w))
    return res
