---
layout: post
date:   2021-12-15
categories: math
author: Guillaume Wang
title: |
  Regularized linear models and the 
  Fenchel-Rockafellar duality theorem (III): 
  Classification with the exponential loss
---

> This is the third of a series of posts on optimization of regularized linear models through the lens of duality.
> See the first one [here](/math/2021/08/27/FRDT_generalities.html).


* This will become a table of contents (this text will be scrapped).
{:toc}

We will continue with the notation from last times, in particular:
- the primal problem is
    $$\label{eq:FRDT_primal} \tag{P}
    \min_{w \in \mathcal{W}} \Psi(w) + \mathcal{L}(V w) =: P(w)
    $$
- the dual problem is
    $$\label{eq:FRDT_dual} \tag{D}
    \max_{a \in \mathcal{Y}^*} - \Psi^*(-V^* a) - \mathcal{L}^*(a) =: D(a).
    $$

In a series of recent works, Ziwei Ji and Matus Telgarsky studied the
optimization of linear models for classification with the exponential
loss (and exponential-like losses), making use of duality arguments
[(Ji, Srebro and Telgarsky, 2021)](http://arxiv.org/abs/2107.00595). Personally I found their
derivations a bit heavy on duality black magic, so I spent a bit of time
to understand what was going on. Unsurprisingly, their notion of dual
variable is exactly the same as the variable $$a$$ of the dual problem
$$\eqref{eq:FRDT_dual}$$. But it actually took me a while to realize
that, for a relatively subtle reason.

In this section, we adopt some of the notation from Ziwei Ji and Matus
Telgarsky's papers, on top of the generic ones already used so far.

-   Let $$\ell(u) = \exp(u)$$ be the exponential loss.

-   We may assume WLOG that $$y^{\text{tgt}}_i = -1$$ for all $$i$$, since
    we may transform the dataset $$(\phi(x_i), y^{\text{tgt}}_i)_i$$ into
    the equivalent dataset $$(z_i, -1)_i$$ with
    $$z_i = -y^{\text{tgt}}_i \phi(x_i)$$. So the data-fitting term
    $$\mathcal{L}(y) = \sum_i \ell(-y^{\text{tgt}}_i y_i)$$, is simply
    $$\mathcal{L}(y) = \sum_i \ell(y_i)$$.

-   The (unnormalized) empirical risk of a parameter $$w$$ is defined as
    $$\mathcal{R}(w) = \sum_{i=1}^n \ell(\left\langle w, z_i \right\rangle)$$.

For classification tasks, it is common to consider as data-fitting term

$$
\mathcal{L}(y) = \sum_{i=1}^n \ell(y_i).
$$

A common trick when
analyzing learning algorithms for classification, is that $$\ell$$ is a
strictly increasing function so that we may use a different choice for
the data-fitting term:
$$\widetilde{\mathcal{L}}= (\ell^{-1}) \circ \mathcal{L}$$, i.e

$$
\widetilde{\mathcal{L}}(y) = \ell^{-1} \left( \sum_{i=1}^n \ell(y_i) \right).
$$

Plus, for $$\ell = \exp$$, $$\widetilde{\mathcal{L}}$$ is just the
log-sum-exp function
$$\widetilde{\mathcal{L}}(y) = \log \sum_{i=1}^n e^{y_i}$$, which is
convex. [^1] [^2]

It turns out that those two seemingly equivalent choices lead to
slightly different optimization algorithms, with significantly different
convergence speeds [(Ji and Telgarsky, 2020)](https://arxiv.org/abs/1906.04540).

1.  Associated with the choice of $$\mathcal{L}$$ is the vanilla gradient
    step 
    
    $$
    w_{t+1} = w_t - \eta_t \nabla \mathcal{R}(w).
    $$

2.  Associated with the choice of $$\widetilde{\mathcal{L}}$$ is the
    normalized gradient step
    
    $$
    w_{t+1} = w_t - \eta_t \frac{\nabla \mathcal{R}(w)}{\mathcal{R}(w)}.
    $$

## Derivation of the two variants of gradient step

Let a convex regularizer $$\Psi(w)$$. Let us naively write down the update
rules from for the saddle-point formulation of the optimization problems
$$\min_w \Psi(w) + \mathcal{L}(Vw)$$ and
$$\min_w \Psi(w) + \widetilde{\mathcal{L}}(Vw)$$. Note that:

-   Since $$\mathcal{L}(Vw) = \mathcal{R}(w)$$,

    $$
    \partial_w \mathcal{L}(Vw) 
            = V^* \partial \mathcal{L}(Vw) 
            = \nabla \mathcal{R}(w).
    $$

-   Since $$\ell^{-1}(v) = \log(v)$$ and so
    $$\widetilde{\mathcal{L}}= \log \circ \mathcal{L}$$,
    
    $$
    \partial_w \widetilde{\mathcal{L}}(Vw)
            = V^* \partial \widetilde{\mathcal{L}}(Vw)
            = \frac{\nabla \mathcal{R}(w)}{\mathcal{R}(w)}.
    $$

Now consider using the scheme from with gradient descent steps for
$$w_{t+1}$$ and fully-optimizing for $$a_{t+1}$$. We get the update rule

$$w_{t+1}
    = w_t - \eta_t \left[ \partial \Psi(w_t) + V^* a_t \right]
    = w_t - \eta_t V^* \partial \mathcal{L}(Vw) - \eta_t \partial \Psi(w_t)
$$

and similarly with $$\widetilde{\mathcal{L}}$$. Plugging in the values of
$$V^* \partial \mathcal{L}(Vw)$$ and
$$V^* \partial \widetilde{\mathcal{L}}(Vw)$$, we see that we get almost
exactly the vanilla and basic gradient steps from above; the only
difference is that we get an extra term $$- \eta_t \partial \Psi(w_t)$$.
When $$\Psi = \frac{\lambda}{2} \left\lVert \cdot \right\rVert_2^2$$, then
as discussed in , a cheap heuristic for implicit regularization (i.e
$$\lambda \to 0$$) is to simply remove that extra term.

It looks like we didn't do anything else than write down the classical
primal gradient descent steps on the unregularized losses
$$\mathcal{L}(Vw)$$ and $$\widetilde{\mathcal{L}}(Vw)$$. That is true. The
advantage of invoking the dual space in this context is that it allows a
finer convergence analysis than if we only stay in the primal
[(Ji and Telgarsky, 2020)](https://arxiv.org/abs/1906.04540). It also leads naturally to a dual accelerated
method, discussed next, that would otherwise appear as utter magic. It
might even allow to divinate yet other funky update rules, by using
other choices for the dual update, or by replacing the exponential by
some other surrogate loss.

#### The dual accelerated method.

In [(Ji, Srebro and Telgarsky, 2021)](http://arxiv.org/abs/2107.00595), they propose a dual-accelerated method for the same
problem (Algorithm 1 of the paper). To present it would require a
discussion of accelerated mirror descent, which would take us a bit far.
Let us only say that their method is essentially just a variant of what
we called the "fully dual approach" [last time](/math/2021/09/03/FRDT_zoo_primal_dual.html#duality-gap-formulation-and-fully-dual-approach-the-frank-wolfe-algorithm), with mirror descent replaced by
a form of accelerated mirror descent.

Interestingly, their new method can also be interpreted as an instance
of the general mix-and-match scheme of , with what seems to be an
unusual form of accelerated gradient descent for $$w_{t+1}$$. However this
point of view is not the one they used to derive and analyze their
method. I find it interesting, and pretty confusing, that a method
derived by acceleration in the dual can be interpreted as a primal-dual
method with acceleration in the primal.

## Beyond $$\ell_2$$: algorithms for $$\ell_1$$-regularized classification, and acceleration

### (F)ISTA for $$\ell_1$$-penalized classification

Consider the same optimization problem as before:
$$\min_w \Psi(w) + \widetilde{\mathcal{L}}(Vw)$$, this time with the
choice of regularizer $$\Psi(w) = \lambda \left\lVert w \right\rVert_1$$.
Consider using the scheme from with fully-optimizing for $$a_{t+1}$$ and
proximal gradient descent for $$w_{t+1}$$. We get the update rule

$$
\begin{aligned}
    a_{t+1} &= \nabla \widetilde{\mathcal{L}}(Vw_t) \\
    w_{t+1}
    &= \mathop{\mathrm{prox}}_{\tau \Psi}(w_t - \tau V^* a_{t+1})\end{aligned}
$$
    
Since $$\Psi = \lambda \left\lVert w \right\rVert_1$$, this is simply the
ISTA algorithm applied to $$\widetilde{\mathcal{L}}(Vw)$$. [^3]

#### Primal acceleration: FISTA.

Accelerated proximal gradient descent in the primal. 

$$\begin{aligned}
    a_{t+1} &= \nabla \widetilde{\mathcal{L}}(V \overline{\gamma}_t) \\
    w_{t+1} &= \mathop{\mathrm{prox}}_{\tau \Psi} \left( \overline{\gamma}_t - \tau V^* a_{t+1} \right) \\
    \overline{\gamma}_{t+1} &= w_{t+1} + \theta (w_{t+1} - w_t)\end{aligned}
$$
    
which reduces to 

$$
\begin{aligned}
    w_{t+1} &= \mathop{\mathrm{prox}}_{\tau \Psi} \left( \overline{\gamma}_t - \tau
    \left.\nabla_w \widetilde{\mathcal{L}}(V w)\right|_{\overline{\gamma}_t}
    \right) \\
    \overline{\gamma}_{t+1} &= w_{t+1} + \theta (w_{t+1} - w_t)\end{aligned}
    
$$

Notice that there is another way to accelerate, by updating $$w_{t+1}$$
starting from $$w_t$$ instead of $$\overline{\gamma}_t$$:

$$\begin{aligned}
    a_{t+1} &= \nabla \widetilde{\mathcal{L}}(V \overline{\gamma}_t) \\
    w_{t+1} &= \mathop{\mathrm{prox}}_{\tau \Psi} \left( w_t - \tau V^* a_{t+1} \right) \\
    \overline{\gamma}_{t+1} &= w_{t+1} + \theta (w_{t+1} - w_t)\end{aligned}
$$

which reduces to 
$$\begin{aligned}
    w_{t+1} &= \mathop{\mathrm{prox}}_{\tau \Psi} \left( w_t - \tau
    \left.\nabla_w \widetilde{\mathcal{L}}(V w)\right|_{\overline{\gamma}_t}
    \right) \\
    \overline{\gamma}_{t+1} &= w_{t+1} + \theta (w_{t+1} - w_t)\end{aligned}
$$
    
This method can be viewed as the Chambolle-Pock algorithm with
$$\sigma = +\infty$$.

### AdaBoost for implicit $$\ell_1$$-regularized classification

It is well-known that AdaBoost results in $$\ell_1$$-margin maximization
[(Chinot, Kuchelmeister, Löffler and van de Geer, 2021)](https://arxiv.org/abs/2105.02083). In this paragraph, we heuristically recover
that fact, by interpreting AdaBoost as (almost) an instance of an
algorithm previously derived in the framework of FRDT.

We consider AdaBoost as stated in Algorithm 1 of
[(Chinot, Kuchelmeister, Löffler and van de Geer, 2021)](https://arxiv.org/abs/2105.02083). [^4] With our notation, one can check that the
algorithm can be formulated as: 

$$\begin{aligned}
    w_0 &= 0 \\
    a_t &= \nabla \widetilde{\mathcal{L}}(V w_t) \\
    u_t &= \left\lVert V^* a_t \right\rVert_\infty \partial \left\lVert \cdot \right\rVert_\infty(-V^* a_t)
    = \partial \frac{1}{2} \left\lVert \cdot \right\rVert_\infty^2 (-V^* a_t) 
    = \partial \left( \frac{1}{2} \left\lVert \cdot \right\rVert_1^2 \right)^* (-V^* a_t) \\
    w_{t+1} &= w_t + \eta u_t\end{aligned}
$$

Denote
$$\varphi: \left[ \mathbb{R}\to \mathbb{R}, x \mapsto \frac{x^2}{2} \right]$$
and $$\psi(w) = \varphi(\left\lVert w \right\rVert_1)$$. In the above
algorithm, the equation for $$u_t$$ can be written as
$$u_t = \partial \psi^*(-V^* a_t)$$, and more generally we have that for
even and convex $$\varphi$$ [^5]

$$\begin{gathered}
    \psi^* = \varphi^* \circ \left\lVert \cdot \right\rVert_\infty, \\
    u_t = \partial \psi^*(-V^* a_t)
    = (\varphi^*)'(\left\lVert -V^* a_t \right\rVert_\infty)~ \partial \left\lVert \cdot \right\rVert_\infty(-V^* a_t).\end{gathered}
$$

So by using different choices for the scalar mapping $$\varphi$$, we
obtain different choices for the adaptive stepsize. We may expect
AdaBoost to have similar regularization behavior for all of them.

Note that AdaBoost is thus strongly reminiscent of the Frank-Wolfe-like
method obtained by what we called the "[fully dual](/math/2021/09/03/FRDT_zoo_primal_dual.html#duality-gap-formulation-and-fully-dual-approach-the-frank-wolfe-algorithm) approach":

$$\begin{aligned}
    w_0, a_0 & ~\text{such that}~ a_0 \in \partial \widetilde{\mathcal{L}}(Vw_0) \\
    a_t &= \nabla \widetilde{\mathcal{L}}(V w_t) \\
    w_{t+1} &= (1-\eta) w_t + \eta \frac{1}{\lambda} \nabla \psi^*(-V^* a_t)\end{aligned}
$$

with $$\lambda \to 0$$. The only difference is that AdaBoost updates
$$w_{t+1}$$ from $$w_t$$ and $$u_t$$ via an additive step instead of a convex
combination. However since $$\lambda \to 0$$, the update
$$\eta \frac{1}{\lambda} \nabla\psi^*(-V^* a_t)$$ can be expected to have
large magnitude so that the $$-\eta w_t$$ term makes no big difference
anyway.

#### Dually accelerated AdaBoost

Armed with the above almost-interpretation of AdaBoost as a previously
derived method, we may derive an accelerated version of AdaBoost. This
would require a discussion of accelerated mirror descent, which would
take us a bit far. Let us only point out that all the necessary
ingredients are contained in Appendix B of [(Ji, Srebro and Telgarsky, 2021)](http://arxiv.org/abs/2107.00595). Namely, I think
the only adaptation needed is to replace $$-V^* a_t$$ ($$-Z^\top q_t$$ in
their notation) by $$u_t = \partial \psi^*(-V^* a_t)$$ everywhere in their
Algorithm 1.

In fact, I expect that deriving and obtaining guarantees for fast
$$\ell_1$$-margin maximization is a very straightforward task, by making
the appropriate adaptations in the proofs of that paper.


[^1]: Beyond the exponential loss, the same trick can be applied for
    other choices of surrogate loss $$\ell$$. A crucial condition for the
    trick is that $$\widetilde{\mathcal{L}}$$ must be convex; additional
    desirable conditions are described in Assumption 1.2 of
    [(Ji and Telgarsky, 2020)](https://arxiv.org/abs/1906.04540), where they also give sufficient
    conditions in their Lemma 5.2.

[^2]: The subtlety that confused me for a while, is that choosing
    $$\mathcal{L}$$ vs. $$\widetilde{\mathcal{L}}$$ as the data-fitting term
    leads to different notions of dual variable,
    $$a = \nabla \mathcal{L}(Vw)$$
    vs. $$\tilde{a}= \nabla \widetilde{\mathcal{L}}(Vw) = (\ell^{-1})'(\mathcal{L}(Vw))~ a$$.
    I initially only had the choice of $$\mathcal{L}$$ in mind, so as I
    stared at the derivations of [(Ji and Telgarsky, 2020)](https://arxiv.org/abs/1906.04540), I could not
    understand why they considered renormalizing the stepsize by
    $$(\ell^{-1})'(\mathcal{L}(Vw_t))$$ in the dual.

[^3]: <https://blogs.princeton.edu/imabandit/2013/04/11/orf523-ista-and-fista/>

[^4]: Our discussion extends immediately to a number of variants of
    AdaBoost: logistic instead of exponential loss, and various choices
    of adaptive stepsize (see the paragraph just below Algorithm 1 in
    [(Chinot, Kuchelmeister, Löffler and van de Geer, 2021)](https://arxiv.org/abs/2105.02083)).

[^5]: See Example 13.8 in the book _Convex Analysis and Monotone Operator Theory in Hilbert Spaces_ by Bauschke and Combettes, 2017.


