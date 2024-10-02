---
layout: post
date:   2021-07-25
last_updated: 2024-09-26
categories: math
author: Guillaume Wang
title: |
  Symmetric tensor functions are dense in the space of
  permutation-invariant functions
---

> I prepared this post a long while ago but only posted it here in July 2021. This is because I had written it in LaTeX and converting it to MD was not completely trivial. Since I had no incentive to post it, this small barrier was enough for me to procrastinate several months...

During my Master's thesis, I encountered several interesting technical
points that were not directly related to the thesis topic, so that I
left them hanging. Here I will talk about one of them, encountered while
working on Volterra series. (I chose to jump directly into my main point
without giving any context on Volterra series, as it is not necessary;
for a clean introduction to these objects, see chapter 4 of my Master's
thesis [report](/contents/Master_s_thesis_report-final.pdf).)


* This will become a table of contents (this text will be scrapped).
{:toc}

### Preliminaries

**Notations and shorthands**
Fix some integer $$n > 0$$.

-   For a point $$\boldsymbol{t}= (t_1,...,t_n) \in \mathbb{R}^n$$ and a
    permutation $$\sigma \in \mathfrak{S}_n$$, $$\boldsymbol{t}_\sigma$$ denotes
    $$(t_{\sigma(1)},...,t_{\sigma(n)})$$.

-   Call a multivariate function $$g: \mathbb{R}^n \to \mathbb{R}$$
    *permutation-invariant* if for any permutation $$\sigma$$, it holds
    $$g(\boldsymbol{t}) = g(\boldsymbol{t}_\sigma)$$ for all $$\boldsymbol{t}\in \mathbb{R}^n$$.

-   For any function $$f: \mathbb{R}\to \mathbb{R}$$, denote
    $$f^{\otimes n}: \left[
            \mathbb{R}^n \to \mathbb{R}, 
            \boldsymbol{t}\mapsto f(t_1)...f(t_n)
        \right]$$. Call
    $$f^{\otimes n}$$ the associated *symmetric tensor function*[^1] --
    tensor because it is a product of single-variable functions, and
    symmetric because all of those single-variable functions are the
    same.

-   For any multivariate function $$g: \mathbb{R}^n \to \mathbb{R}$$,
    denote $$\mathop{\mathrm{Sym}}g: \left[
            \mathbb{R}^n \to \mathbb{R},
            \boldsymbol{t}\mapsto \frac{1}{n!} \sum_{\sigma \in \mathfrak{S}_n} g(\boldsymbol{t}_\sigma)
        \right]$$.

We will sometimes write physicist-style $$f(t)$$ to mean a function
$$f: \mathbb{R}\to \mathbb{R}$$, and similarly $$g(\boldsymbol{t})$$ instead of
$$g: \mathbb{R}^n \to \mathbb{R}$$.

**Some function spaces**
Fix $$1 \leq p < \infty$$ and $$q$$ its conjugate exponent, i.e., $$1/p+1/q=1$$.

-   Let $$L^p(\mathbb{R})$$ be the Banach space of $$L^p$$-integrable
    functions over $$\mathbb{R}$$ (with the usual Lebesgue measure). Its
    dual space is $$L^q(\mathbb{R})$$.

-   Let $$C_0(\mathbb{R})$$ be the space of vanishing continuous functions
    over $$\mathbb{R}$$. Its dual space is $$\mathcal{M}(\mathbb{R})$$, the
    space of Radon measures.[^2]

-   Similarly define $$L^p(\mathbb{R}^n)$$, $$L^q(\mathbb{R}^n)$$,
    $$C_0(\mathbb{R}^n)$$ spaces of multivariate functions with $$n$$ scalar
    variables.

-   Denote $$L^p_{\mathop{\mathrm{Sym}}}(\mathbb{R}^n)$$,
    $$L^q_{\mathop{\mathrm{Sym}}}(\mathbb{R}^n)$$,
    $$C_{0 \mathop{\mathrm{Sym}}}(\mathbb{R}^n)$$ the respective (closed)
    subspaces consisting of permutation-invariant functions.

Note that our shorthand $$\mathop{\mathrm{Sym}}$$ can be viewed as a
projection operator from $$L^p(\mathbb{R}^n)$$ to
$$L^p_{\mathop{\mathrm{Sym}}}(\mathbb{R}^n)$$, and from
$$C_0(\mathbb{R}^n)$$ to $$C_{0 \mathop{\mathrm{Sym}}}(\mathbb{R}^n)$$.

### The result and why it looks surprising to me

<div class="proposition" markdown="1" text="main result">
Let $$1 \leq p < \infty$$. The set $$\left\lbrace
        f^{\otimes n}(\boldsymbol{t}) ;~ f \in L^p(\mathbb{R})
    \right\rbrace$$ has its linear span dense in
$$L^p_{\mathop{\mathrm{Sym}}}(\mathbb{R}^n)$$.

The set $$\left\lbrace
        f^{\otimes n}(\boldsymbol{t}) ;~ f \in C_0(\mathbb{R})
    \right\rbrace$$ has its linear span dense in
$$C_{0 \mathop{\mathrm{Sym}}}(\mathbb{R}^n)$$.[^3]
</div>

More explicitly: *any
$$g(\boldsymbol{t}) \in C_{0 \mathop{\mathrm{Sym}}}(\mathbb{R}^n)$$ is
arbitrarily-well uniformly approximated by finite sums of the form
$$\sum_{i \leq m} f_i(t_1)...f_i(t_n)$$* ($$m < \infty$$,
$$f_i \in C_0(\mathbb{R})$$).

**This is not Weierstrass with symmetrization**
As an obvious corollary, the proposition holds when $$\mathbb{R}$$ is
replaced by a closed interval $$I \subset \mathbb{R}$$. In this case the
result looks like a straightforward consequence of the Weierstrass
approximation theorem, but it is not. Consider the following valid
reasoning:

> Fix a continuous function $$g(\boldsymbol{t})$$ over the compact $$I^n$$ and let
> $$\varepsilon>0$$. By the Weierstrass approximation theorem, there
> exists a polynomial $$P(\boldsymbol{t})$$ such that
> $$\left\lVert g-P \right\rVert := \sup_{I^n} \left\lvert g-P \right\rvert \leq \varepsilon$$,
> and $$P(\boldsymbol{t})$$ can be written as
> $$P(\boldsymbol{t}) = \sum_{\alpha \in \mathbb{N}^n} a_\alpha \boldsymbol{t}^\alpha$$
> (where there are only a finite number of nonzero coefficients
> $$a_\alpha$$ and the shorthand $$\boldsymbol{t}^\alpha$$ denotes
> $$t_1^{\alpha_1} ... t_n^{\alpha_n}$$).
>
> If in addition $$g(\boldsymbol{t})$$ is permutation-invariant, then the lemma
> below shows that $$\mathop{\mathrm{Sym}}P(\boldsymbol{t})$$ is also an
> $$\varepsilon$$-approximation of $$g(\boldsymbol{t})$$, and it can be written as
> $$
> \mathop{\mathrm{Sym}}P(\boldsymbol{t}) = \sum_{\alpha \in \mathbb{N}^n} a_\alpha \mathop{\mathrm{Sym}}\boldsymbol{t}^\alpha
>         = \sum_{\alpha \in \mathbb{N}^n} \sum_{\sigma \in \mathfrak{S}_n} \frac{a_\alpha}{n!} t_1^{\alpha_{\sigma(1)}} ... t_n^{\alpha_{\sigma(n)}}.
> $$

This does not show that $$g(\boldsymbol{t})$$ can be approximated by a finite
combination of symmetric tensor functions, as the reasoning may yield
approximators such as $$t_1 t_2^3 + t_1^3 t_2$$ (if $$n=2$$), which are not
of the required form.

<div class="lemma" markdown="1" text='the aforementioned lemma'>
For any function $$g(\boldsymbol{t})$$ over $$\mathbb{R}^n$$ and any
$$1 \leq p \leq \infty$$, it holds:
$$\left\lVert \mathop{\mathrm{Sym}}g \right\rVert_{L^p} \leq \left\lVert g \right\rVert_{L^p}$$.

For any permutation-invariant $$g(\boldsymbol{t})$$ and any function $$h(\boldsymbol{t})$$,
if $$\left\lVert g - h \right\rVert_{L^p} \leq \varepsilon$$, then
$$\left\lVert g - \mathop{\mathrm{Sym}}h \right\rVert_{L^p} \leq \varepsilon$$.
</div>

<div class="proof" markdown="1">
Let any function $$g(\boldsymbol{t})$$ over $$\mathbb{R}^n$$ and any
$$1 \leq p \leq \infty$$. By definition of the $$L^p$$ norm,

$$
\left\lVert \mathop{\mathrm{Sym}}g \right\rVert_{L^p} = \left\lVert  \frac{1}{n!} \sum_{\sigma \in \mathfrak{S}_n} g(\boldsymbol{t}_\sigma)  \right\rVert_{L^p}
        \leq \frac{1}{n!} \sum_{\sigma \in \mathfrak{S}_n} \left\lVert g(\boldsymbol{t}_\sigma) \right\rVert_{L^p}
        = \left\lVert g \right\rVert_{L^p}.
$$


Let any permutation-invariant function $$g(\boldsymbol{t})$$ and any function
$$h(\boldsymbol{t})$$ such that
$$\left\lVert g - h \right\rVert_{L^p} \leq \varepsilon$$. Then

$$
\left\lVert g - \mathop{\mathrm{Sym}}h \right\rVert_{L^p} = \left\lVert \mathop{\mathrm{Sym}}(g-h) \right\rVert_{L^p} \leq \left\lVert g-h \right\rVert_{L^p} \leq \varepsilon.
$$
</div>


**An example of surprise**
In fact the case of permutation-invariant polynomials over a compact set
is already surprising to me\... In fact just the following example is
already surprising to me:

> Consider the function $$g(\boldsymbol{t}) = t_1 + ... + t_n$$ over $$[0,1]^n$$.
> According to the proposition, there exist $$m<\infty$$ and
> $$f_i(t) \in C([0,1])$$ ($$i \leq m$$) such that
> $$g(\boldsymbol{t}) \approx \sum_{i \leq m} f_i(t_1)...f_i(t_n)$$, in the sense
> of uniform approximation over $$[0,1]^n$$.

I wonder what these $$f_i$$ could look like. Since they are continuous
over $$[0,1]$$, according to the Weierstrass approximation theorem we may
assume without loss of generality that each $$f_i$$ is polynomial. Then,
developing the product $$f_i(t_1)...f_i(t_n)$$ would yield an a priori big
polynomial, whereas directly using Weierstrass with symmetrization can
yield simply $$t_1+...+t_n$$ itself.

### Brief proof of the result

In this section we prove the $$L^p/L^q$$ ($$1 \leq p < \infty$$) part of the
proposition; the $$C_0/\mathcal{M}$$ part can be proved by the same
arguments, with minor modifications.

I assume the reader is familiar with the basics of functional analysis
and duality in Banach spaces. Recall the following density criterion,
which is a consequence of the Hahn-Banach theorem (as are many things):

<div class="lemma" markdown="1" text="density criterion">
Let $$E$$ be a Banach space and $$A$$ a subset. 
$$
\begin{aligned}
        \mathop{\mathrm{span}}(A) \text{ is dense in } E
        && \iff &&
        \left\lbrace
            X \in E';~
            \forall a \in A, \left\langle a, X \right\rangle_E = 0
        \right\rbrace
        = \{ 0_{E'} \}.
    \end{aligned}
$$
</div>

The set on the right is sometimes denoted $$A^\perp$$ and called the
*annihilator* of $$A$$. Note that
$$A^\perp = (\mathop{\mathrm{span}}(A))^\perp$$. In the case where $$E$$ is
a Euclidean space, $$A^\perp$$ is just (up to isometry) the orthogonal
complement of $$\mathop{\mathrm{span}}(A)$$.

The proposition will be proved by applying the above density criterion.
To do so we will need the following intuitively obvious lemma,
characterizing the dual of $$L^p_{\mathop{\mathrm{Sym}}}(\mathbb{R}^n)$$.
A formal and rather uninteresting proof can be found
in appendix of the [LaTeX version](/assets/pdf/symmetric_funs_density.pdf).

<div class="lemma" markdown="1">
The dual of $$L^p_{\mathop{\mathrm{Sym}}}(\mathbb{R}^n)$$ is isometrically
isomorphic to $$L^q_{\mathop{\mathrm{Sym}}}(\mathbb{R}^n)$$.
</div>

We can now prove the proposition. The main argument is extracted from
(Boyd Chua Desoer 1984, Theorem 2.5.2) in the context of Volterra
series.

<div class="proof" markdown="1" text="of proposition">
To apply the density criterion to $$A = \left\lbrace
        f^{\otimes n}(\boldsymbol{t}) ;~ f \in L^p(\mathbb{R})
    \right\rbrace$$ and $$E = L^p_{\mathop{\mathrm{Sym}}}(\mathbb{R}^n)$$,
let $$h \in E' \simeq L^q_{\mathop{\mathrm{Sym}}}(\mathbb{R}^n)$$ such
that $$\left\langle f^{\otimes n}, h \right\rangle_{L^p} = 0$$ for all
$$f \in L^p(\mathbb{R})$$. Let us show that $$h=0$$, from which the
proposition will follow.

Denote $$\Phi_h: L^p(\mathbb{R}) \to \mathbb{R}$$ the $$n$$-homogeneous map

$$
\Phi_h[f]
        = \left\langle f^{\otimes n}, h \right\rangle_{L^p}
        = \int_{\mathbb{R}^n} d\boldsymbol{t}~ h(t_1,...,t_n) f(t_1)...f(t_n)
$$

and $$\Psi_h: L^p(\mathbb{R})^n \to \mathbb{R}$$ the associated $$n$$-linear
system 

$$
\Psi_h\{f_1,...,f_n\}
        = \left\langle f_1 \otimes ... \otimes f_n, h \right\rangle_{L^p}
        = \int_{\mathbb{R}^n} d\boldsymbol{t}~ h(t_1,...,t_n) f_1(t_1)...f_n(t_n).
$$


The $$n$$-linear system $$\Psi_h\{\cdot,...,\cdot\}$$ is symmetric in its
arguments, since $$h$$ is permutation-invariant. So $$\Psi_h$$ is completely
determined by the $$n$$-homogeneous map $$\Phi_h[\cdot]$$ via the [algebraic
polarization identity](https://en.wikipedia.org/wiki/Polarization_of_an_algebraic_form)

$$
n! \Psi_h\{f_1,...,f_n\} 
        = \left.
            \frac{\partial}{\partial \alpha_1 ... \partial \alpha_n}
        \right|_{\boldsymbol{\alpha}=0}
        \Phi_h \left[ \sum_{i=1}^n \alpha_i f_i \right],
$$

and the
right-hand-side is the differential of an identically zero map.
Consequently,

$$
\forall f_1,...,f_n \in L^p(\mathbb{R}),~ \Psi_h\{f_1,...,f_n\} = 0.
$$


Now evaluate this at
$$f_1(t) = \boldsymbol{1}_{t \in A_1}, ..., f_n(t) = \boldsymbol{1}_{t \in A_n}$$
for intervals $$A_i \subset \mathbb{R}$$:

$$
\Psi_h\{f_1,...,f_n\} = \int_{\mathbb{R}^n} d\boldsymbol{t}~ h(t_1,...,t_n)
        \boldsymbol{1}_{\boldsymbol{t}\in A_1 \times ... \times A_n} = 0.
$$

Since this holds for all $$A_i$$, and hyperrectangles generate the Borel
$$\sigma$$-algebra, then $$h=0$$, as claimed.
</div>

### Is the result interesting/useful?

For the subjects that I'm currently leaning towards, the result
presented in this document is actually pretty useless, as it only talks
about approximability per se. It doesn't give any guarantees on the
nature nor the number of functions $$f_i$$ required to
$$\varepsilon$$-approximate a given target function $$g$$.

However I still find the result technically interesting and surprising.
I never heard about it before but I'm certain it must be somewhere out
there already -- I would be glad to know where and in what context.

------

[^1]: Disclaimer: the term \"symmetric tensor function\" may not be
    consistent with standard terminology, I haven't checked.

[^2]: <https://regularize.wordpress.com/2011/11/11/dual-spaces-of-continuous-functions/>

[^3]: I'm pretty sure the same holds if $$C_0$$ is replaced by $$C_b$$, i.e,
    if we consider bounded continuous functions, instead of vanishing
    continuous.
