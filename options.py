def shapeOf(x):
    """
    By convention if x is an array, the shape should be n*dim
    """
    if type(x) in [int, float]:
        return 1, 1
    else:
        X = np.array(x)
        try :
            n, dim = X.shape
        except:
            n = 1
            dim = X.shape[0]

    return n, dim