#ifndef LAB_NONLINEARSOLVE_H
#define LAB_NONLINEARSOLVE_H

std::vector<TT> Newton(const std::string &method, size_t dim, const std::vector<TT> &y0);

std::vector<TT> f(const std::vector<TT>& x);

#endif // LAB_NONLINEARSOLVE_H
