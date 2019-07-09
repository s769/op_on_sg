from plotting import plot_monomial, plot_SOP

# 1. To plot the monomials, use plot_monomial
# num represent the highest power to be plotted
# k represent the type of monomial {0,1,2} -> {1, 2, 3}
num = 5
k = 0

plot_monomial(num, k)

# 2. This plots the Sobolov Orthogonal Polynomials
# num represent the highest power to be plotted
# k represent the type of monomial {0,1,2} -> {1, 2, 3}
plot_SOP(num, k)