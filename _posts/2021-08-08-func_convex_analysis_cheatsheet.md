---
layout: post
date:   2021-08-08
categories: math
author: Guillaume Wang
title: A functional and convex analysis cheat sheet
---

> Another post that was prepared a while ago and "snoozed" until now...
> It is part of a planned series of posts on linear models and regularization, and a tinge of optimization.

Among the topics that I've been interested in during the last few
months, many required some knowledge of convex analysis: regularization
in linear models, primal-dual views on optimization, representer
theorems in reproducing kernel Banach spaces\... Moreover the
finite-dimensional setting doesn't suffice, a more abstract point of
view is necessary or at least useful; namely Banach spaces seem to be
the appropriate level of abstraction for those topics.

In this document I compile some relevant functional and convex analysis
background, in the form of a cheat sheet. It is not at all meant to be
exhaustive, I only included basic facts and tricks that I found
interesting. I may add to it in the future.

Proofs and appendices can be found in the [LaTeX version](/contents/func_convex_analysis_cheatsheet.pdf) of this document.

* This will become a table of contents (this text will be scrapped).
{:toc}

## Functional analysis (Banach duality)

Beyond finite dimension, Banach spaces are a simple and natural level of
abstraction for discussing convex analysis. This section is mostly
extracted from the appendix of my [Master's thesis](/contents/Master_s_thesis_report-final.pdf).

<div class="definition" markdown="1" text="Banach space">
A metric space $$(E,d)$$ is called *complete* if all Cauchy sequences
$$(u_n)_n \in E^\mathbb{N}$$ converge in $$E$$.

A *Banach space* $$(E, \left\lVert \cdot \right\rVert_E)$$ is a vector
space equipped with a norm for which it is a complete space.

A *Hilbert space* $$(H, \left\langle \cdot, \cdot \right\rangle_H)$$ is a
vector space equipped with an inner product that is complete for the
induced norm
$$\left\lVert x \right\rVert_H^2 = \left\langle x, x \right\rangle_H$$.

The unit ball of a normed space $$(E, \left\lVert \cdot \right\rVert_E)$$
is denoted $$B^{(E)} = B^{(E)}_{0,1} 
    := \left\lbrace x \in E;~ \left\lVert x \right\rVert_E \leq 1 \right\rbrace$$.

A continuous linear mapping $$T: E \to F$$ between Banach spaces is called
a *bounded operator*, and its *operator norm* is the finite quantity
$${\left\vert\kern-0.25ex\left\vert\kern-0.25ex\left\vert T \right\vert\kern-0.25ex\right\vert\kern-0.25ex\right\vert} = {\left\vert\kern-0.25ex\left\vert\kern-0.25ex\left\vert T \right\vert\kern-0.25ex\right\vert\kern-0.25ex\right\vert}_{E \to F}
    := \sup_{\left\lVert x \right\rVert_E \leq 1} \left\lVert T x \right\rVert_F$$.
The set of bounded operators from $$E$$ to $$F$$ equipped with the operator
norm
$$(\mathcal{L}_b(E,F), {\left\vert\kern-0.25ex\left\vert\kern-0.25ex\left\vert \cdot \right\vert\kern-0.25ex\right\vert\kern-0.25ex\right\vert})$$
is itself a Banach space.

A bounded operator $$T: E \to F$$ is called *compact* if it sends the unit
ball into a relatively compact set, i.e $$T(B^{(E)})$$ is a relatively
compact set of $$F$$, i.e $$\overline{T(B^{(E)}})$$ is compact where
$$\overline{~\cdot~}$$ denotes closure w.r.t the norm of $$F$$.
</div>

### Duality in Banach spaces

<div class="definition" markdown="1" text="dual space">
The *(topological) dual* of a Banach space $$E$$ is the space of bounded
linear forms $$E' = \mathcal{L}_b(E, \mathbb{R})$$. It is equipped with
the norm
$$\left\lVert X \right\rVert_{E'} := \sup_{\left\lVert x \right\rVert_E \leq 1} \left\lvert X(x) \right\rvert$$.
$$E'$$ is itself a Banach space.

The *duality bracket* of $$E$$ is the bilinear operator
$$\left\langle \cdot, \cdot \right\rangle_E: E \times E' \to \mathbb{R}$$
defined by $$\left\langle x, X \right\rangle_E = X(x)$$.

The *bidual* of $$E$$ is the space $$E'' = (E')'$$. $$E$$ can be embedded into
$$E''$$ by $$x \mapsto x''$$, where $$x''$$ is defined by:
$$\forall X \in E',~ \left\langle X, x'' \right\rangle_{E'} = \left\langle x, X \right\rangle_E = X(x)$$.
It is not hard to show (using existence of norming functionals, see
below) that this embedding is isometric i.e
$$\left\lVert x \right\rVert_{(E')'} = \left\lVert x \right\rVert_E$$.

$$E$$ is called a *reflexive* Banach space if the converse holds, i.e if
any element of the bidual can also be seen as an element of the primal,
i.e if $$E'' \simeq E$$.
</div>

In this section, elements of the primal space will typically be denoted
by lowercase letters e.g $$x \in E, y \in F$$, and elements of the dual by
uppercase letters e.g $$X \in E', Y \in F'$$.

<div class="remark" markdown="1" text="bra-ket">
The duality bracket is very similar to the physicists' bra-ket notation;
except that here the primal is on the left and the dual is on the right,
instead of the opposite.

When $$E$$ is reflexive, then all the shorthands from the bra-ket notation
can be used. That is, a dual element $$X$$ can be denoted without
ambiguity as $$\left\langle \cdot, X \right\rangle_E$$ , and a primal
element $$x = x''$$ as $$\left\langle x, \cdot \right\rangle_E$$ . Moreover,
for a bounded operator $$T: E \to F$$, we can write without ambiguity
$$\left\langle Tx, Y \right\rangle_F = \left\langle x| T |Y \right\rangle$$ .
However since there are many interesting Banach spaces that are not
reflexive, we will not use such shorthands.
</div>

### Hahn-Banach theorem and useful consequences

<div class="theorem" markdown="1" text="Hahn-Banach">
Let $$\underline{E}$$ be a linear subspace of a normed vector space
$$(E, \left\lVert \cdot \right\rVert)$$ and let
$$f: \underline{E}\to \mathbb{R}$$ be a bounded linear form on
$$(\underline{E}, \left\lVert \cdot \right\rVert)$$.

Then there exists $$g: E \to \mathbb{R}$$ a bounded linear form on all of
$$E$$, such that

-   $$g$$ is an extension of $$f$$:
    $$\left.g\right|_{\underline{E}} = f$$;

-   The extension "comes to no cost" in operator norm:
    $${\left\vert\kern-0.25ex\left\vert\kern-0.25ex\left\vert g \right\vert\kern-0.25ex\right\vert\kern-0.25ex\right\vert} = {\left\vert\kern-0.25ex\left\vert\kern-0.25ex\left\vert f \right\vert\kern-0.25ex\right\vert\kern-0.25ex\right\vert}$$.

(Here the operator norms are with respect to their respective domains:
$${\left\vert\kern-0.25ex\left\vert\kern-0.25ex\left\vert f \right\vert\kern-0.25ex\right\vert\kern-0.25ex\right\vert} = \sup_{x \in \underline{E}; \left\lVert x \right\rVert \leq 1} \left\lvert f(x) \right\rvert$$,
$${\left\vert\kern-0.25ex\left\vert\kern-0.25ex\left\vert g \right\vert\kern-0.25ex\right\vert\kern-0.25ex\right\vert} = \sup_{x \in E; \left\lVert x \right\rVert \leq 1} \left\lvert g(x) \right\rvert$$.)
</div>

As one of the many important consequences of that theorem, we have the
existence of norming functionals.

<div class="definition" markdown="1" text="norming functional">
Let $$E$$ be a Banach space.

For all $$x \in E \setminus \{0_E\}$$, there exists $$X \in E'$$ such that
$$X(x) = \left\langle x, X \right\rangle_E = \left\lVert x \right\rVert_E$$
and $$\left\lVert X \right\rVert_{E'} = 1$$. $$X$$ is then called a *norming
functional* of $$x$$.

By convention, any $$X \in B^{(E')}$$ will be called a norming functional
of $$0_E$$.
</div>

Importantly, the norming functional $$X$$ is not unique in general, and
there is no generic way to construct it -- the proof is not
constructive. [^1] (This stands in contrast with the case of Hilbert
spaces, where the norming functional is unique and given by the Riesz
representation theorem.)

<div class="proposition" markdown="1">
The primal $$E$$ injects isometrically into the bidual $$E''$$.
Namely, $$x \in E$$ is canonically associated to $$x'' \in E''$$ defined by
$$\forall X \in E', \left\langle X, x'' \right\rangle_{E'} = \left\langle x, X \right\rangle_E$$.
</div>

Another important consequence of the Hahn-Banach theorem is the
following density criterion.

<div class="lemma" markdown="1" text="Hahn-Banach density criterion">
Let $$E$$ be a Banach space. Let $$\underline{E}$$ be any subspace and $$A$$
any subset. 

$$
\begin{aligned}
    \underline{E}\text{ is dense in } E
    && \iff &&
    \underline{E}^\perp :=
    \left\lbrace
        X \in E';~
        \forall x \in \underline{E}, \left\langle x, X \right\rangle_E = 0
    \right\rbrace
    = \{ 0_{E'} \} \\
    \mathop{\mathrm{span}}(A) \text{ is dense in } E
    && \iff &&
    A^\perp :=
    \left\lbrace
        X \in E';~
        \forall a \in A, \left\langle a, X \right\rangle_E = 0
    \right\rbrace
    = \{ 0_{E'} \}
\end{aligned}
$$
</div>

The set $$A^\perp$$ is called the *annihilator* of $$A$$. Note that
$$A^\perp = (\mathop{\mathrm{span}}(A))^\perp$$. In the case where $$E$$ is
a Euclidean space, $$A^\perp$$ is just (up to isometry) the orthogonal
complement of $$\mathop{\mathrm{span}}(A)$$.

## Convex analysis (convex duality)

For the rest of this subsection, fix a Banach space $$E$$.

<div class="definition" markdown="1">
A function $$f: E \to \mathbb{R}\cup \{+\infty\}$$ is called *convex* if

$$
\forall x, y \in E, \forall t \in [0,1],~ f \left( tx + (1-t) y \right) \leq t f(x) + (1-t) f(y).
$$

For any convex $$f: E \to \mathbb{R}\cup \{+\infty\}$$,

-   The *domain* of $$f$$ is the convex set
    $$\mathop{\mathrm{dom}}(f) = \left\lbrace 
                x \in E; f(x) < \infty
            \right\rbrace$$. $$f$$ is called *proper* if
    $$\mathop{\mathrm{dom}}(f) \neq \varnothing$$.

-   $$f$$ is called *lower-semicontinuous* (l.s.c) if its sub-level sets
    are closed, i.e for each $$c \in \mathbb{R}$$,
    $$\left\lbrace x \in E; f(x) > c \right\rbrace$$ is an open set.

Denote $$\Gamma(E)$$ the set of proper l.s.c convex functions over $$E$$.
</div>

<div class="definition" markdown="1">
For any proper convex function $$f: E \to \mathbb{R}\cup \{+\infty\}$$,

-   The *subdifferential* of $$f$$ at a point $$x_0 \in E$$ is the set

    $$\partial f(x_0) = \left\lbrace
                    X \in E'; \forall x \in E, f(x) \geq f(x_0) + \left\langle x-x_0, X \right\rangle_E
                \right\rbrace.$$ 
                
    $$f$$ is called *subdifferentiable* at
    $$x_0$$ if $$\partial f(x_0) \neq \varnothing$$. Note that $$f$$ can be
    subdifferentiable at $$x_0$$ only if
    $$x_0 \in \mathop{\mathrm{dom}}(f)$$.

-   $$f$$ is called *differentiable* at $$x_0$$ if $$\partial f(x_0)$$ is a
    singleton. Its unique element is then called the *differential* of
    $$f$$ at $$x_0$$ and denoted $$D f(x_0)$$ or $$\nabla f(x_0)$$.
</div>

Note that the definitions of "subdifferential" and "differential" above
look different from the usual ones, since they only apply to convex
functions. It can be shown that our definitions are compatible with the
usual ones from real analysis. [^2]

<div class="proposition" markdown="1">
For any proper l.s.c function $$f: E \to \mathbb{R}\cup \{+\infty\}$$ such
that $$\mathop{\mathrm{dom}}(f)$$ is open,

-   $$f$$ is convex iff $$\mathop{\mathrm{dom}}(f)$$ is convex and
    $$\partial f(x_0) \neq \varnothing$$ for all
    $$x_0 \in \mathop{\mathrm{dom}}(f)$$.

-   If $$f$$ is convex, then it is differentiable at $$x_0$$ (in the usual
    real-analytic sense) iff $$\partial f(x_0)$$ is a singleton, and the
    differential of $$f$$ at $$x_0$$ (in the usual real-analytic sense) is
    then $$D f(x_0)$$.
</div>

### Convex conjugate

<div class="definition" markdown="1">
For any proper function $$f: E \to \mathbb{R}\cup \{+\infty\}$$, the
*convex conjugate* of $$f$$ (a.k.a Fenchel-Legendre a.k.a Legendre-Fenchel
a.k.a Fenchel a.k.a Legendre transform) is the function

$$
f^*: \left[ E' \to \mathbb{R}\cup \{+\infty\},
        X \mapsto \sup_{x \in E} \left\langle x, X \right\rangle_E - f(x) \right].
$$
</div>

<div class="proposition" markdown="1" text="Fenchel-Moreau theorem">
For any proper function $$f$$, $$f^*$$ is a proper l.s.c convex function.

A function $$f$$ is a proper l.s.c convex function iff $$f^{**} = f$$.
</div>

For any proper function $$f$$, $$f^{**}$$ is the tightest convex relaxation
of $$f$$, in the sense that the epigraph of $$f^{**}$$ is the convex hull of
the epigraph of $$f$$. This is easy to visualize for functions over the
real line.

<div class="remark" markdown="1">
In the proposition above, $$f^{**}$$ is understood as a mapping from $$E$$
to $$\mathbb{R}$$. Looking at the definitions, it would be more natural to
view $$f^{**}$$ as a mapping from $$E''$$ to $$\mathbb{R}$$ instead, which
would be more general since $$E$$ injects isometrically into $$E''$$.
However in convex analysis we typically don't care about what happens
outside of $$E$$.

More precisely: to be completely general and consistent with notation,
we could define $$f^{**} = (f^*)^*$$ over $$E''$$ by
$$
\forall z \in E'',~ f^{**}(z) = \sup_{X \in E'} \left\langle X, z \right\rangle_{E'} - f^*(X).
$$

Since
$$\left\langle X, x'' \right\rangle_{E'} = \left\langle x, X \right\rangle_E$$
(where $$\left[ E \to E'', x \mapsto x'' \right]$$ denotes the canonical
injection), the restriction of $$f^{**}$$ to $$E$$ is then -- and this
equation is typically taken as the definition of $$f^{**}$$:

$$
\forall x \in E,~ f^{**}(x) =
        \sup_{X \in E'} \left\langle x, X \right\rangle_E - f^*(X).
$$
</div>

In the context of convex analysis it is common to denote
adjoint/conjugate/dual objects with a superscript \"$$*$$\". In other
contexts that symbol connotes involution, which may be misleading.
However for convex analysis there is not much risk of mistake, precisely
because of the previous remark: we only ever care about what happens in
$$E$$ and $$E'$$, never about the bidual space $$E''$$. In particular, even if
$$f$$ is not a proper l.s.c function, we may always write
$$(f^{**})^* = f^*$$.

Accordingly, **from here on we will follow the common practice and use
$$x^*$$ (instead of $$X$$) to denote a generic element of $$E'$$.**

Many of the useful properties of convex conjugates can be found on the relevant [wikipedia page](https://en.wikipedia.org/w/index.php?title=Convex_conjugate&oldid=1007941296), so I won't list those here again.

### Convex conjugates vs. subdifferentials

<div class="proposition" markdown="1" text="Fenchel-Young inequality">
Let $$f \in \Gamma(E)$$. By definition of the convex conjugate, we have
Fenchel-Young's inequality: 

$$\forall x \in E, \forall x^* \in E',~
        f(x) + f^*(x^*) \geq \left\langle x, x^* \right\rangle_E.$$

For
$$x \in E$$ and $$x^* \in E'$$,

-   $$x^* \in \partial f(x)$$ iff $$x^*$$ saturates Fenchel-Young's
    inequality, iff $$x^*$$ achieves the sup in the definition of
    $$f(x) = f^{**}(x) = \sup_{x^* \in E'} \left\langle x, x^* \right\rangle_E - f^*(x^*)$$.

-   $$x \in \partial f^*(x^*)$$ iff $$x$$ saturates Fenchel-Young's
    inequality, iff $$x$$ achieves the sup in the definition of
    $$f^*(x^*) = \sup_{x \in E} \left\langle x, x^* \right\rangle_E - f(x)$$.

-   $$x^* \in \partial f(x)$$ iff $$x \in \partial f^*(x^*)$$.
</div>

<div class="remark" markdown="1" text="subdifferentials as correspondence, Rockafellar 1970, Theorem 24.9">
Up to an additive constant,
$$f \in \Gamma(E)$$ is characterized by the binary relation $$\mathcal{R}$$
given by

$$
x \mathcal{R}x^* \iff f(x) + f^*(x^*) = \left\langle x, x^* \right\rangle_E.
$$

</div>

<div class="proposition" markdown="1" text="norming functionals as subdifferentials">
The norm $$\left\lVert \cdot \right\rVert_E$$ is a proper continuous
convex function by definition.

For any $$x \in E$$, $$x^* \in E'$$ is a norming functional for $$x$$ iff
$$x^* \in \partial \left\lVert \cdot \right\rVert_E(x)$$. In symbols,

$$
x^* \in \partial \left\lVert \cdot \right\rVert_E(x)
\iff
\begin{cases}
    \left\lVert x^* \right\rVert_{E'} = 1 \\
    \left\langle x, x^* \right\rangle = \left\lVert x \right\rVert_E
\end{cases}
$$

In particular, $$\left\lVert \cdot \right\rVert_E$$ is differentiable at
$$x$$ iff $$\partial \left\lVert \cdot \right\rVert_E(x)$$ is a singleton,
iff $$x$$ has a unique norming functional. If
$$\left\lVert \cdot \right\rVert_E$$ is differentiable everywhere, then
the mapping
$$\left[ x \mapsto \nabla \left\lVert \cdot \right\rVert_E(x) \right]$$ is
well-defined and is called the *duality mapping*.
</div>

Thus, to the convex analyst, norming functionals are not a magical
byproduct of the Hahn-Banach theorem, but simply a subgradient of the
norm.

### Convex conjugacy swaps strict convexity for differentiability, and strong convexity for smoothness

<div class="definition" markdown="1">
Let $$f \in \Gamma(E)$$ and let $$\mu>0$$, $$L>0$$.

$$f$$ is *strictly convex* if for all $$x_0 \in \mathop{\mathrm{dom}}(f)$$,
there exists $$g \in \partial f(x_0)$$ such that the strict inequalities
hold:

$$\forall x \in E \setminus \{x_0\},~ f(x) > f(x_0) + \left\langle x-x_0, g \right\rangle_E.$$

$$f$$ is [*$$\mu$$-strongly convex*](https://xingyuzhou.org/blog/notes/strong-convexity) if for all
$$x_0 \in \mathop{\mathrm{dom}}(f)$$, there exists $$g \in \partial f(x_0)$$
such that

$$\forall x \in E,~ f(x) \geq f(x_0) + \left\langle x-x_0, g \right\rangle_E + \frac{\mu}{2} \left\lVert x-x_0 \right\rVert_E^2.$$

$$f$$ is differentiable everywhere if it is differentiable at each point
of its domain, that is, if for all $$x_0 \in \mathop{\mathrm{dom}}(f)$$,
there exists *a unique* $$g \in E'$$ such that

$$\forall x \in E,~ f(x) \leq f(x_0) + \left\langle x-x_0, g \right\rangle_E.$$

$$f$$ is *$$L$$-smooth* if for all $$x_0 \in \mathop{\mathrm{dom}}(f)$$, there
exists $$g \in \partial f(x_0)$$ such that

$$\forall x \in E,~ f(x) \leq f(x_0) + \left\langle x-x_0, g \right\rangle_E + \frac{L}{2} \left\lVert x-x_0 \right\rVert_E^2.$$

Note that $$\mu$$-strong convexity implies strict convexity, and that
$$L$$-smoothness implies differentiability everywhere.
</div>

<div class="proposition" markdown="1" text="Kakade, Shalev-Shwartz, and Tewari, 2009">
Let $$f \in \Gamma(E)$$.

$$f$$ is strictly convex iff $$f^*$$ is differentiable.

$$f$$ is $$\mu$$-strongly convex iff $$f^*$$ is $$1/\mu$$-smooth.
</div>

------

[^1]: Or rather, to the pure functional analyst there is no satisfying
    way to construct a norming functional, but to the convex analyst
    there is\... See below.

[^2]: In other words, in this document we are only concerned with
    content covered in Rockafellar's 1970 "Convex Analysis" book,
    whereas in other contexts Rockafellar's 1998 "Variational Analysis"
    book may be a better reference.
