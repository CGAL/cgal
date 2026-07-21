# LBFGS-Lite
A header-only L-BFGS unconstrained optimizer.

## 0. About

__LBFGS-Lite__ is a __C++__ [__header-only__](https://en.wikipedia.org/wiki/Header-only) library for __unconstrained optimization__. Many engineering considerations are added for improved robustness compared to the original version by Nocedal.

## 1. How to use

All explanations are detailed by the comments in "lbfgs.hpp". See "lbfgs_example.cpp" for the calling procedure. You may need to install Eigen via "sudo apt install libeigen3-dev" because we use Eigen for better performance since [ver. 2.1](https://github.com/ZJU-FAST-Lab/LBFGS-Lite/tags). If you need a pure C-style lib without any external dependence, please refer to [ver. 0.9](https://github.com/ZJU-FAST-Lab/LBFGS-Lite/tags).

## 2. Features

- Only one header file "lbfgs.hpp" is all you need.

- The library implements [Limited-Memory Broyden-Fletcher-Goldfarb-Shanno Method](https://doi.org/10.1007/BF01589116) (L-BFGS).

- A __highly robust line search__ proposed by [Lewis and Overton](https://link.springer.com/article/10.1007/s10107-012-0514-2) has been employed since [ver. 2.1](https://github.com/ZJU-FAST-Lab/LBFGS-Lite/tags).

- Both __smooth__ ([C2](https://en.wikipedia.org/wiki/Smoothness)) and __nonsmooth__ ([C0 but piecewise C2](https://en.wikipedia.org/wiki/Smoothness)) functions are supported.

- __Cautious update__ by [Li and Fukushima](https://epubs.siam.org/doi/pdf/10.1137/S1052623499354242) is employed for __global convergence__ in nonconvex cases.

- __Externally provided maximum step size__ is convenient to functions defined on bounded sets.

## 3. Why this lib

- Our lib is well-organized and concise (about 800 lines).

- Many other libs use [Nocedal's zoom line search](https://link.springer.com/book/10.1007%2F978-0-387-40065-5), [More-Thuente line search](https://dl.acm.org/doi/abs/10.1145/192115.192132), or [Hager-Zhang line search](https://doi.org/10.1137/030601880). The truth is, interpolation-based schemes bring many tunable parameters and make more or less assumptions on the function smoothness. Engineering practice tells us that these assumptions can fail due to bad numerical conditions and ill-shaped functions. Admittedly, they slightly reduce the iteration number in some cases, but also easily-crashed.

- Some other libs provide user-specified options for [Armijo and weak/strong Wolfe condition](https://en.wikipedia.org/wiki/Wolfe_conditions) without considering positive definiteness (PD) of the approximated Hessian. This is misleading because if only Armijo condition is satisfied, the PD property is not guaranteed and the solver is unstable or easily-crashed. We make the weak Wolfe condition mandatory here, which suffices for PD property.

- We use [Lewis-Overton line search](https://link.springer.com/article/10.1007/s10107-012-0514-2) as the only scheme since ver. 2.1 from which nonsmooth functions are supported. Other schemes either assume high orders of continuity, or enforce the strong Wolfe condition can never be fulfilled by nonsmooth functions. Moreover, Lewis-Overton line search are widely adopted in many graphics applications involving optimization on [scene](https://dl.acm.org/doi/abs/10.1145/2766929), [shape](https://dl.acm.org/doi/abs/10.1145/2897824.2925918), or [mesh](https://dl.acm.org/doi/abs/10.1145/3197517.3201303), showing its practical robustness.

- According to our practice, the function/gradient evaluation numbers are comparable with interpolation-based schemes. Sometimes it even requires fewer function calls. If you insist an interpolation-one for smooth well-shaped cost function, we also propose our [ver. 0.9](https://github.com/ZJU-FAST-Lab/LBFGS-Lite/tags) where a More-Thuente line search is kept.

- Other schemes' global convergence on non-convex functions are not guaranteed theoretically. We avoid the potential problems by employing the [cautious update](https://epubs.siam.org/doi/pdf/10.1137/S1052623499354242) scheme in our lib without additional computation overhead.

## 6. Licence

LBFGS-Lite is modified from [the C version](https://github.com/chokkan/liblbfgs) by Okazaki, which is further based on [the original Fortran version](https://doi.org/10.1007/BF01589116) by Nocedal. Thus it is distributed under the term of the MIT license according to previous ones by Okazaki and by Nocedal. Please refer to LICENSE file in the distribution.

## 7. Maintaince

If any bug, please contact [Zhepei Wang](https://zhepeiwang.github.io/) (<wangzhepei@live.com>).
