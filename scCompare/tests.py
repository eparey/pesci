






def init_ec(mat1, mat2):
    ec = []
    for i in range(len(mat1[:,0])): #avoid looping in later optimization
        ec.append(max(np.corrcoef(mat1[i,:], mat2[i,:])[0,1], 0))
    return np.array(ec)



def init_ec_optimized(mat1, mat2, n):
    
    #row-wise correlation and retain argmax



    for i in range(len(mat1[:,0])): #avoid looping in later optimization
        ec.append(max(np.corrcoef(mat1[i,:], mat2[i,:])[0,1], 0))
    return np.array(ec)




#covariance / (variance1 * variance2)
#   n = len(mat1[:,0])
def corrcoef_einsum(mat1, mat2, n):

    avg1 = np.einsum("ij->i", mat1, optimize='optimal') / np.double(n)
    avg2 = np.einsum("ij->i", mat2, optimize='optimal') / np.double(n)

    mat1_m = mat1 - avg1
    mat2_m  = mat2 - avg2

    var1 = np.einsum("ij,ij->j", mat1_m, mat1_m, optimize='optimal')
    var2 = np.einsum("ij,ij->j", mat2_m, mat2_m, optimize='optimal')

    cov = np.einsum("ij,ik->jk", mat1_m, mat2_m, optimize='optimal')
    tmp = np.einsum("i,j->ij", var1, var2, optimize='optimal')
    res = np.diag(cov / np.sqrt(tmp), 0).clip(min=0)
    return res


def corrcoef_einsum2(mat1, mat2, n):

    mat1_m = mat1 - np.einsum("ij->j", mat1, optimize='optimal') / np.double(n)
    mat2_m = mat2 - np.einsum("ij->j", mat2, optimize='optimal') / np.double(n)

    cov = np.einsum("ij,ij->j", mat1_m, mat2_m, optimize='optimal')

    var1 = np.einsum("ij,ij->j", mat1_m, mat1_m, optimize='optimal')
    var2 = np.einsum("ij,ij->j", mat2_m, mat2_m, optimize='optimal')
    varprod = np.einsum("i,i->i", var1, var2, optimize='optimal')

    res = cov / np.sqrt(varprod)
    return res.clip(min=0)


def wcov_einsum(mat1, mat2, weights, n):
 
    wmat = np.tile(weights, (n, 1))
    mat1_m = mat1 - np.einsum("ji,ij->j", wmat, mat1, optimize='optimal') / sum(weights)
    mat2_m = mat2 - np.einsum("ji,ij->j", wmat, mat2, optimize='optimal') / sum(weights)
 
    cov = np.einsum("ij,ij->ij", mat1_m, mat2_m, optimize='optimal')
    wcov = np.einsum("ji,ij->j", wmat, cov, optimize='optimal') / sum(weights)
 
    cov1 = np.einsum("ij,ij->ij", mat1_m, mat1_m, optimize='optimal')
    wcov1 = np.einsum("ji,ij->j", wmat, cov1, optimize='optimal') / sum(weights)
 
    cov2 = np.einsum("ij,ij->ij", mat2_m, mat2_m, optimize='optimal')
    wcov2 = np.einsum("ji,ij->j", wmat, cov2, optimize='optimal') / sum(weights)
 
    prod = np.einsum("i,i->i", wcov1, wcov2, optimize='optimal')
 
    res = wcov / np.sqrt(prod)
    return res.clip(min=0)



def wcov(v1, v2, w):
    return np.sum(w * (v1 - np.average(v1, weights=w)) * (v2 - np.average(v2, weights=w))) / np.sum(w)
    #np.average((v1 - np.average(v1, weights=w)) * (v2 - np.average(v2, weights=w)), weights=w)

def wcorr(v1, v2, w):
    return wcov(v1, v2, w) / np.sqrt(wcov(v1, v1, w) * wcov(v2, v2, w))