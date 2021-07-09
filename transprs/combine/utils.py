from scipy.optimize import nnls


def nonneg_lstsq(X, y):

    assert np.all(X.std(axis=0) > 0)
    y_mean = y.mean()
    X_mean = X.mean(axis=0)
    y_c = y - y_mean
    X_c = X - X_mean
    coef, _ = nnls(X_c, y_c)
    y_hat = X.dot(coef)
    intercept = np.mean(y - y_hat)
    return coef, intercept
