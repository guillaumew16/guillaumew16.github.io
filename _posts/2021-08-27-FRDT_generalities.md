---
layout: post
date:   2021-08-27
categories: math
author: Guillaume Wang
title: |
  Regularized linear models and the 
  Fenchel-Rockafellar duality theorem (I):
  Generalities
---

> This is the first of a series of posts on optimization of regularized linear models through the lens of duality.
> The series will be essentially a summary of what I learned during the first few weeks of my internship, when I looked into topics related to learning in Banach spaces (before I moved on to different, more concrete topics).
> The relevant convex analysis background can be found in last time's [cheatsheet](/blog/2021/func_convex_analysis_cheatsheet) -- which is basically the zero-th post of this series.
> 
> With the background definitions under our belt, we are *almost* ready to talk about concrete consequences of convex duality.
> As it turns out, to get a principled understanding of the many places where duality pops up, it is beneficial to first present the Fenchel-Rockafellar duality theorem, a kind of "master theorem" for convex duality.


Many smart observations on (explicitly or implicitly) regularized linear
models, as well as on convex optimization methods, make use of some sort
of convex duality argument. For example, already in introductory ML
courses, Lagrangian duality is typically used to show that the SVM
solution is equivalent to the max-margin classifier. In this document we
present a unified framework for these duality arguments, that makes it
conceptually easier to draw connections.

On a personal note, I have always found convex and Langrangian duality
to be particularly black magic. The nice thing about the
Fenchel-Rockafellar duality theorem is that it is general enough to
contain all of that black magic, making duality-based derivations easier
to follow.

* This will become a table of contents (this text will be scrapped).
{:toc}

### Generic supervised learning setup and notation {#sec:setup_notation}

Consider supervised learning with a linear model (linear in the
parameters) i.e a hypothesis space of the form

$$
\mathcal{F}= \left\lbrace
    f_w: x \mapsto \left\langle w, \phi(x) \right\rangle,
    w \in \mathcal{W}
\right\rbrace
$$

for some feature map
$$\phi: \mathcal{X}\to \mathcal{W}^*$$, where $$\mathcal{X}$$ is input space
and $$\mathcal{W}$$ is a Banach parameter space.

Suppose we are given a dataset $$(x_i, y^{\text{tgt}}_i)_{i \leq n}$$ and
we want to solve an optimization problem of the form [^1]

$$
\mathop{\mathrm{arg\,min}}_{w \in \mathcal{W}} \Psi(w) + \frac{1}{n} \sum_{i=1}^n \ell \left( y^{\text{tgt}}_i, \left\langle w, \phi(x_i) \right\rangle \right),
$$

where each $$\ell(y^{\text{tgt}}_i, \cdot)$$ is a convex function, and the
regularizer term $$\Psi$$ is convex. This generic formulation captures
pretty much all standard supervised learning settings, see the last section of this post.

Further pose the shorthands:

-   $$\mathcal{Y}= (\mathbb{R}^n, \left\lVert \cdot \right\rVert)$$ for
    some arbitrary norm; $$\mathcal{Y}^*$$ will denote $$\mathbb{R}^n$$
    equipped with the dual norm $$\left\lVert \cdot \right\rVert_*$$;

-   $$\mathcal{L}(y) = \frac{1}{n} \sum_i \ell(y^{\text{tgt}}_i, y_i)$$
    the data-fitting term;

-   $$V: \left[ \mathcal{W}\to \mathcal{Y}, w \mapsto (\left\langle w, \phi(x_i) \right\rangle)_{i \leq n} \right]$$
    the evaluation operator.

#### Standard notations for convex duality.

-   For a Banach space $$E$$, the dual space is denoted $$E^*$$.

-   The set of proper, lower-semicontinuous (l.s.c), and convex
    functions over $$E$$ is denoted $$\Gamma(E)$$.

-   For a convex function $$\Psi$$, $$\Psi^*$$ denotes the convex conjugate.

-   For a linear operator $$V: \mathcal{W}\to \mathcal{Y}$$ between Banach
    spaces, $$V^*: \mathcal{Y}^* \to \mathcal{W}^*$$ denotes the adjoint
    operator between the dual spaces.

-   The indicator function of a set $$A$$ is defined as
    $$\iota_A(x) = \begin{cases}
            \infty \text{ if } x \not\in A \\
            0 \text{ if } x \in A
        \end{cases}$$.

-   $$\left\lVert \cdot \right\rVert$$ denotes an arbitrary norm, and
    $$\left\lVert \cdot \right\rVert_*$$ denotes its dual norm.

### Fenchel-Rockafellar duality theorem (FRDT) {#sec:FRDT}

Let us state the FRDT with machine-learning-friendly notation, as
motivated above.

<div class="theorem" markdown="1" text="FRDT">
Let $$\mathcal{W}$$ and $$\mathcal{Y}$$ be two real Banach spaces. Let
$$\Psi \in \Gamma(\mathcal{W})$$, let
$$\mathcal{L}\in \Gamma(\mathcal{Y})$$, and let
$$V: \mathcal{W}\to \mathcal{Y}$$ be a bounded linear operator. Consider
the primal problem 

$$\label{eq:FRDT_primal} \tag{P}
\min_{w \in \mathcal{W}} \Psi(w) + \mathcal{L}(V w) =: P(w)
$$

and define the dual problem as 

$$\label{eq:FRDT_dual} \tag{D}
\max_{a \in \mathcal{Y}^*} - \Psi^*(-V^* a) - \mathcal{L}^*(a) =: D(a).
$$

Denote $$S_P$$ and $$S_D$$ their respective optimal solution sets. Denote
the KKT conditions 

$$\label{eq:FRDT_KKT_Psi} \tag{$\mathrm{KKT}_\Psi$}
w \in \partial \Psi^*(-V^* a)
~~\text{i.e}~~
-V^* a \in \partial \Psi(w)
$$

$$\label{eq:FRDT_KKT_L} \tag{$\mathrm{KKT}_{\mathcal{L}}$}
a \in \partial \mathcal{L}(V w)
~~\text{i.e}~~
V w \in \partial \mathcal{L}^*(a)
$$

Weak duality holds: for all $$w$$ and $$a$$,
$$P(w) \geq P_{\text{opt}} \geq D_{\text{opt}} \geq D(a)$$.

Suppose that
$$0 \in \text{interior}( V(\mathop{\mathrm{dom}}\Psi)- \mathop{\mathrm{dom}}\mathcal{L})$$.
Then strong duality holds:

-   $$\eqref{eq:FRDT_primal}$$ and $$\eqref{eq:FRDT_dual}$$ have the same optimal value
    $$P_{\text{opt}} = D_{\text{opt}}$$;

-   $$w \in S_P$$ and $$a \in S_D$$ iff
    $$\eqref{eq:FRDT_KKT_L}$$ and $$\eqref{eq:FRDT_KKT_Psi}$$;

-   $$w \in S_P$$ iff there exists $$a \in \mathcal{Y}^*$$ such that
    $$\eqref{eq:FRDT_KKT_Psi}$$ and
    $$\eqref{eq:FRDT_KKT_L}$$, iff there exists $$a \in S_D$$ such that
    $$\eqref{eq:FRDT_KKT_L}$$.
</div>

<div class="proof" markdown="1" text="sketched">
The actual proof is much more complicated, [^2] but hinges on the
following calculation which is enough to gain intuition:

$$
\begin{aligned}
    \eqref{eq:FRDT_primal}
    &\equiv \min_w \Psi(w) + \mathcal{L}(Vw) \\
    &\equiv \min_w \max_a \Psi(w) + \left\langle Vw, a \right\rangle_\mathcal{Y}- \mathcal{L}^*(a) \\
    &\equiv \min_w \max_a \Psi(w) + \left\langle w, V^* a \right\rangle_\mathcal{W}- \mathcal{L}^*(a) \\
    &\geq \max_a \min_w - \left[ \left\langle w, -V^* a \right\rangle_\mathcal{Y}- \Psi(w) \right] - \mathcal{L}^*(a) \\
    &\equiv -\Psi^*(-V^*a) - \mathcal{L}^*(a)
    \equiv \eqref{eq:FRDT_dual}
\end{aligned}
$$

Denote
$$F(w,a) = \Psi(w) + \left\langle Vw, a \right\rangle - \mathcal{L}^*(a)$$.
The above calculation shows that $$P(w) \geq F(w,a) \geq D(a)$$ for all
$$w,a$$. To see where the KKT conditions come from, note that
$$P(w) = D(a)$$ implies

-   $$F(w,a) = P(w)$$ i.e
    $$\left\langle Vw, a \right\rangle_\mathcal{Y}- \mathcal{L}^*(a)= \mathcal{L}(Vw)$$,
    i.e $$a$$ saturates the Fenchel-Young inequality, i.e
    $$a \in \partial \mathcal{L}(Vw)$$;

-   $$F(w,a) = D(a)$$ i.e
    $$\left\langle w, -V^* a \right\rangle_\mathcal{Y}- \Psi(w) = \Psi^*(-V^*a)$$,
    i.e $$w$$ saturates the Fenchel-Young inequality, i.e
    $$w \in \partial \Psi^*(-V^* a)$$.
</div>

<div class="remark" markdown="1">
The condition that
$$0 \in \text{interior}( V(\mathop{\mathrm{dom}}\Psi)- \mathop{\mathrm{dom}}\mathcal{L})$$
is morally just a variant of Slater's condition: "there exists a
strictly feasible point". There are other constraint qualification
conditions that imply strong duality.
</div>

<div class="remark" markdown="1">
The adjoint of the evaluation operator, $$V^*$$, is simply given by

$$
\forall a \in \mathcal{Y}^* = \mathbb{R}^n,~
V^* a = \sum_{i=1}^n a_i \left\langle \cdot, \phi(x_i) \right\rangle
= \left\langle \cdot, \sum_{i=1}^n a_i \phi(x_i) \right\rangle.
$$

For finite-dimensional features, say $$\dim(\mathcal{W}) = p$$, $$V$$ can be
seen as the transformed data matrix $$\begin{bmatrix}
        \phi(x_1) & ... & \phi(x_n)
    \end{bmatrix}^\top \in \mathbb{R}^{n \times p}$$, and $$V^*$$ is simply
its transpose.
</div>

### Variants of regularization: penalized, constrained, min-norm interpolation {#sec:variants_of_regu}

It is common wisdom that *models with lower "complexity" have better
generalization properties*.

Traditionally, low "complexity" of the learned model is ensured in
practice by adding a penalty term to the loss function. In our notation:
$$f_w$$ is chosen as minimizing $$\lambda \psi(w) + \mathcal{L}(Vw)$$, for
some regularizer or "complexity measure" $$\psi$$. Still traditionally, a
theory-friendlier alternative is to constrain the parameters to a
low-complexity set $$\Omega$$, e.g
$$\Omega = \left\lbrace w; \psi(w) \leq B \right\rbrace$$.

Parallel to these ideas, also of interest for theory and practice are
so-called *overparametrized* settings, whereby
$$\dim(\mathcal{W}) \gg \dim(\mathcal{X})$$, so that there exist many
values of $$w$$ that interpolate the data i.e such that
$$\mathcal{L}(V w) = 0$$. In such settings, the aforementioned common
wisdom can be interpreted in a third way: select, among all
interpolating values of $$w$$, the one with the lowest "complexity":
$$\mathop{\mathrm{arg\,min}}_w \psi(w) ~\text{s.t}~ \mathcal{L}(V w) = 0$$.
[^3]

Note that, since FRDT strictly generalizes Lagrangian duality, it offers
a unified framework for:

-   regularization via penalization, by letting $$\Psi(w)$$ be the penalty
    term $$\lambda \psi(w)$$;

-   regularization via constraining the parameters to a convex set
    $$\Omega$$, by letting $$\Psi(w) = \iota_\Omega(w)$$;

-   min-regularizer e.g min-norm interpolation, by letting
    $$\mathcal{L}(Vw) = \iota_{ \{y^{\text{tgt}}\} }(Vw)$$.

Regardless of the setting (penalized, constrained, or interpolation),
many smart observations on linear models can be made by using some sort
of duality argument. The nice thing about the FRDT is that it provides a
unified way to formulate all those duality-based arguments.

#### Coming up next...

In the next few posts, I will discuss a few example topics where
convex duality arguments are invoked, through the lens of the FRDT.
I hope to convince you that adopting that viewpoint is indeed very helpful
for getting a principled understanding of those topics.
Up next: a zoo of primal-dual methods for optimization of regularized linear models.

------

[^1]: Usually the dataset is denoted $$(x_i, y_i)_{i \leq n}$$ and the
    predictions are denoted $$\hat{y}_i$$. Here we chose to denote
    $$y^{\text{tgt}}_i$$ the target labels because it will be more
    convenient to denote simply $$y_i$$ the predictions.

[^2]: The proof can be found in Section 31 of Rockafellar's 1970 book
    "Convex analysis", and on this webpage:
    <https://pwacker.com/fenchelrockafellar.html>, which also contains
    nice illustrative drawings. Apparently this blog post series was
    also planning to nicely present the proof:
    <https://dohmatob.github.io/research/2019/10/31/duality.html>, but
    it hasn't been updated in a while.

[^3]: Note that this formalization of \"data-interpolation\" is only for
    regression, since for classification with the logistic loss for
    example, $$\mathcal{L}(V w) = 0$$ is impossible.
    