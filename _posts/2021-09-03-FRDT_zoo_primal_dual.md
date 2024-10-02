---
layout: post
date:   2021-09-03
categories: math
author: Guillaume Wang
title: |
  Regularized linear models and the 
  Fenchel-Rockafellar duality theorem (II):
  A zoo of primal-dual methods
---

> This is the second of a series of posts on optimization of regularized linear models through the lens of duality.
> See the first one [here](/blog/2021/FRDT_generalities).


* This will become a table of contents (this text will be scrapped).
{:toc}

We will continue with the notation from last time, in particular:
- the primal problem is
    
    $$\label{eq:FRDT_primal} \tag{P}
    \min_{w \in \mathcal{W}} \Psi(w) + \mathcal{L}(V w) =: P(w)
    $$

- the dual problem is
    
    $$\label{eq:FRDT_dual} \tag{D}
    \max_{a \in \mathcal{Y}^*} - \Psi^*(-V^* a) - \mathcal{L}^*(a) =: D(a).
    $$

- the KKT conditions are
    
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


## Saddle-point formulation: "mix-and-match" primal and dual updates

The FRDT essentially tells us that the problem of fitting a regularized
linear model to data, the minimization problem
$$\eqref{eq:FRDT_primal}$$, can be formulated as a saddle-point
problem:

$$\min_{w \in \mathcal{W}} \max_{a \in \mathcal{Y}^*}~
    F(w,a) =:
    \Psi(w) + \left\langle Vw, a \right\rangle - \mathcal{L}^*(a)$$

Let us naively think about how to iteratively solve that saddle-point
problem; that is, how to choose an update rule for the joint variable
$$(w_t, a_t)$$. I can think of four reasonable update rules for $$w_{t+1}$$,
given $$(w_t,a_t)$$:

-   **Fully optimize for fixed $$a$$:** By definition, for a fixed value
    of $$a$$, the $$\mathop{\mathrm{arg\,min}}$$ of the objective over $$w$$
    is computable in closed form as

    $$\mathop{\mathrm{arg\,min}}_{w'} F(w',a) 
            = \mathop{\mathrm{arg\,min}}_{w'} \Psi(w') + \left\langle w', V^* a \right\rangle
            = \mathop{\mathrm{arg\,max}}_{w'} \left\langle w', -V^* a \right\rangle  - \Psi(w')
            = \partial \Psi^*(-V^* a).
    $$
    
    So we may take as update rule:
    
    $$w_{t+1} \in \partial \Psi^*(-V^* a_t)$$
    
    (This can be interpreted as enforcing the condition
    $$\eqref{eq:FRDT_KKT_Psi}$$ throughout the optimization procedure.)
    However, this update rule may not be computationally feasible.

-   **Gradient descent step:** Instead of fully optimizing, we may take
    one (or several) (sub)gradient descent step(s) for
    $$\mathop{\mathrm{arg\,min}}_{w'} F(w',a_t)$$ starting from $$w_t$$.
    This gives the update rule:

    $$w_{t+1} \in w_t - \eta_t \partial_w F(w_t,a_t)
            = w_t - \eta_t \left[ \partial \Psi(w_t) + V^* a_t \right]$$

-   **Mirror descent step using $$\Psi$$:** Instead of gradient descent,
    we may take one step of the other basic optimization primitive that
    is mirror descent. A natural candidate for the link function is
    $$\Psi$$ (assuming it is strictly convex and differentiable
    everywhere), yielding the update rule

    $$\nabla \Psi(w_{t+1}) = \nabla \Psi(w_t) - \eta_t \left[ \nabla \Psi(w_t) + V^* a_t \right]$$

-   **Proximal gradient step:** Instead of mirror descent, we may take
    one step of the third basic optimization primitive that is proximal
    gradient descent. It is arguably more natural than mirror descent
    since the objective if composite. This gives the update rule

    $$w_{t+1} = \mathrm{prox}_{\eta_t \Psi} \left(
                w_t - \eta_t V^* a_t
            \right)
    $$

As for update rules for $$a_{t+1}$$ given $$(w_t,a_t)$$, the four same ideas
apply.

-   **Fully optimize for fixed $$w$$:**
    $$a_{t+1} \in \partial \mathcal{L}(V w_t)$$

-   **Gradient descent step:**
    $$a_{t+1} \in a_t + \sigma_t \left[ - \partial \mathcal{L}^*(a_t) + V w_t \right]$$

-   **Mirror descent step using $$\mathcal{L}^*$$:**
    $$\nabla \mathcal{L}^*(a_{t+1}) = \nabla \mathcal{L}^*(a_t) + \sigma_t \left[ - \nabla \mathcal{L}^*(a_t) + V w_t \right]$$

-   **Proximal gradient step:**
    $$a_{t+1} = \mathrm{prox}_{\sigma_t \mathcal{L}^*} \left(
                a_t + \sigma_t V w_t
            \right)$$

Thus, this naive reasoning gives a set of optimization schemes that one
can try: simply mix-and-match the choice of update rules for $$w_{t+1}$$
and for $$a_{t+1}$$, and apply the updates alternatingly. By alternating
updates we mean: $$w_{t+1} = \text{Update}_w(w_t,a_t)$$,
$$a_{t+1} = \text{Update}_a(w_{t+1},a_t)$$. One can even consider using
joint updates: $$w_{t+1} = \text{Update}_w(w_t,a_t)$$,
$$a_{t+1} = \text{Update}_a(w_t,a_t)$$.


### Examples

- **Fully-optimizing in the dual recovers vanilla primal methods.**
    Indeed if we substitute $$a_t$$ by $$\partial \mathcal{L}(V w_t)$$ in the
    primal update rules, then the term $$V^* a_t$$ becomes

    $$V^* \partial \mathcal{L}(V w_t) 
        = \left.
            \partial_w \mathcal{L}(V w)
        \right|_{w_t}$$
        
    the subgradient of the data-fitting term w.r.t the primal variable.

- **Proximal gradient steps in the dual.**
    The algorithm consisting of alternating proximal gradient steps both for
    $$w_{t+1}$$ and for $$a_{t+1}$$, is called the Arrow-Hurwicz method.
    The well-known Chambolle-Pock algorithm
    can be interpreted as a fancier version of this scheme, whereby proximal
    gradient steps for $$a_{t+1}$$ are alternated with a form of accelerated
    proximal gradient for $$w_{t+1}$$: [(Chambolle and Pock, 2011)](https://hal.archives-ouvertes.fr/hal-00490826/document)

    $$\begin{aligned}
        a_{t+1} &= \mathrm{prox}_{\sigma \mathcal{L}^*} \left( a_t + \sigma V {\overline{w}}_t \right) \\
        w_{t+1} &= \mathrm{prox}_{\tau \Psi} \left( w_t - \tau V^* a_{t+1} \right) \\
        {\overline{w}}_{t+1} &= w_{t+1} + \theta (w_{t+1} - w_t)\end{aligned}$$

    and the parameters $$\sigma, \tau, \theta$$ can further be made to depend
    on $$t$$.

    As a second example, the proximal dual coordinate ascent algorithm
    proposed in [(Raj and Bach, 2021)](https://arxiv.org/abs/2003.13807), and its accelerated variant, are also
    instances of the scheme explained above. There, the primal variables are
    fully optimized ($$w_t = \nabla \Psi^*(-V^* a_t)$$), and the dual
    variables are updated by a proximal gradient step. That paper focuses on
    the specific case of min-$$\Psi$$-interpolation, so $$\mathcal{L}^*$$
    consists in a sum of indicator functions, and they use an explicit
    expression for $$\mathrm{prox}_{\mathcal{L}^*}$$.

- **Kernel methods (RKHS).**
    For (Hilbert) kernel methods, $$\mathcal{W}$$ is the RKHS and the
    regularizer is
    $$\Psi(w) = \frac{\lambda}{2} \left\lVert w \right\rVert^2$$. So
    $$\partial \Psi^*(-V^* a_t) = -\lambda V^* a_t$$, and one may
    fully-optimize in the primal and run e.g gradient descent entirely in
    the dual. In function space this corresponds to parametrizing the model
    as $$f = \sum_{i=1}^n a_i k(\cdot, x_i)$$ and running gradient descent on
    the coefficients $$a_i$$.



## Duality-gap formulation and fully dual approach: the Frank-Wolfe algorithm

For a given pair of variables $$(w,a)$$, we call *duality gap* the
quantity $$P(w) - D(a)$$. Since
$$P(w) \geq P_{\text{opt}} \geq D_{\text{opt}} \geq D(a)$$, the duality
gap is non-negative, and provides an optimality certificate for the
primal: $$P(w) - P_{\text{opt}} \leq P(W) - D(a)$$. So instead of solving
the saddle-point problem $$\min_w \max_a F(w,a)$$, we may consider solving
the duality-gap minimization problem

$$\min_{w \in \mathcal{W}} \min_{a \in \mathcal{Y}^*} P(w) - D(a).$$

Observe that the duality gap can be split into two terms as:

$$\begin{aligned}
    P(w) - D(a)
    &= \left[ P(w) - F(w,a) \right]
    + \left[ F(w,a) - D(a) \right] \\
    &= \left[ \mathcal{L}^*(a) + \mathcal{L}(Vw) - \left\langle Vw, a \right\rangle \right]
    + \left[ \Psi^*(-V^*a) + \Psi(w) - \left\langle w, -V^* a \right\rangle \right].\end{aligned}$$

Both bracketed terms are non-negative. The first term is zero iff
$$\eqref{eq:FRDT_KKT_L}$$ is satisfied, and the second term is zero iff
$$\eqref{eq:FRDT_KKT_Psi}$$ is satisfied.

The above suggests the following idea:

-   Choose an update rule for $$w_{t+1}$$ such that
    $$\eqref{eq:FRDT_KKT_L}$$ is always satisfied, so that at each
    step, $$F(w_t,a_t) = P(w_t)$$;

-   Choose an update rule for $$a_{t+1}$$ that takes a step towards
    minimizing the second term in the duality gap:
    $$\min_{a'} F(w_t,a') - D(a') = \Psi^*(-V^*a') + \Psi(w_t) - \left\langle w_t, -V^* a' \right\rangle$$.

Note that, compared to the saddle-point paradigm from the previous
subsection, this idea seems completely backwards:

-   We saw that fully-optimizing in the primal for the saddle-point
    problem $$\min_w \max_a F(w,a)$$ leads to choosing $$w_t$$ that always
    satisfies the KKT condition for $$\Psi$$
    $$\eqref{eq:FRDT_KKT_Psi}$$; whereas here we enforce the KKT
    condition for $$\mathcal{L}$$
    $$\eqref{eq:FRDT_KKT_L}$$.

-   For the dual update, the saddle-point approach suggests to choose
    $$a_{t+1}$$ as taking a step towards $$\max_{a'} F(w_t,a')$$, i.e to
    take a step towards minimizing the first term in the duality gap:
    $$\min_{a'} P(w_t) - F(w_t,a')$$; whereas here we take a step towards
    minimizing the second term.

This is due to the fact that here we choose $$w_t$$ to optimize only the
first term of the duality-gap split: $$P(w) - F(w,a)$$, and boldly ignored
the second term...

### A trick to enforce the "wrong" KKT condition

To actually implement the idea explained above, we are faced with a
difficulty. The primal update rule is to choose $$w_t$$ such that
$$V w_t \in \partial \mathcal{L}^*(a_t)$$. A naive approach is to
dumbly compute $$\partial \mathcal{L}^*(a_t)$$ and to somehow find a $$w_t$$
in its preimage by $$V$$. But $$V$$, the evaluation operator, is typically
difficult to invert (think of $$V$$ as the data matrix and
$$V \in \mathbb{R}^{n \times p}$$ with $$p \gg n$$).

Now comes a trick: suppose we have, at timestep $$t$$,
$$V w_t \in \partial \mathcal{L}^*(a_t)$$, and we want to construct
$$w_{t+1}$$ such that $$V w_{t+1} \in \partial \mathcal{L}^*(a_{t+1})$$.
Further suppose that we chose to update $$a_{t+1}$$ by mirror descent
using $$\mathcal{L}^*$$: 

$$\begin{aligned}
    \nabla \mathcal{L}^*(a_{t+1})
    &= \nabla \mathcal{L}^*(a_t) - \sigma_t 
    \left.
    \partial_a \left( 
        \Psi^*(-V^*a) + \Psi(w_t) - \left\langle w_t, -V^* a \right\rangle
    \right)
    \right|_{a_t} \\
    &= \nabla \mathcal{L}^*(a_t) 
    + \sigma_t V \nabla \Psi^*(-V^* a_t)
    - \sigma_t V w_t.\end{aligned}$$
    
Thus, we wish to construct $$w_{t+1}$$ such that 

$$V w_{t+1} = V w_t 
    + \sigma_t V \nabla \Psi^*(-V^* a)
    - \sigma_t V w_t.$$
    
Clearly a possible choice is to simply set
$$w_{t+1} = (1-\sigma_t) w_t + \sigma_t \nabla \Psi^*(-V^* a_t)$$ !

Thus, we obtain the following algorithm for implicit regularization:

$$\begin{aligned}
    w_0, a_0 & ~\text{such that}~ a_0 \in \partial \mathcal{L}(Vw_0) \\
    \nabla \mathcal{L}^*(a_{t+1}) &= \nabla \mathcal{L}^*(a_t) 
    + \sigma_t V \nabla \Psi^*(-V^* a_t)
    - \sigma_t V w_t \\
    w_{t+1} &= (1-\sigma_t) w_t + \sigma_t \nabla \Psi^*(-V^* a_t)\end{aligned}$$

which can be simplified as 

$$\begin{aligned}
    w_0, a_0 & ~\text{such that}~ a_0 \in \partial \mathcal{L}(Vw_0) \\
    a_t &= \nabla \mathcal{L}(V w_t) \\
    w_{t+1} &= (1-\sigma_t) w_t + \sigma_t \nabla \Psi^*(-V^* a_t)\end{aligned}$$

To avoid confusion, note that the first equation (defining $$a_{t+1}$$) is
what we called the primal update rule; and that the second equation
(which looks like a gradient step for $$w_{t+1}$$) is actually the
mirror descent update for the dual.

This trick that allows to choose $$w_{t+1}$$ satisfying
$$\eqref{eq:FRDT_KKT_L}$$, can also be applied to other choices of the
dual update. A crucial ingredient is that the dual update should involve
mirror descent using $$\mathcal{L}^*$$. In particular, the trick can be
applied for accelerated mirror descent in the dual: see the paragraph
just below Lemma 3.2 in [(Ji, Srebro and Telgarsky, 2021)](http://arxiv.org/abs/2107.00595), as well as their Appendix B.

### Relation to Frank-Wolfe

Note that the above method seems significantly different from anything
we could arrive to by mixing-and-matching primal and dual updates for
the saddle-point formulation. Indeed, the primal update rule

$$\begin{aligned}
    u_t &= \partial \Psi^*(-V^*a_t) \\
    w_{t+1} &= (1-\sigma_t) w_t + \sigma_t u_t\end{aligned}$$
    
looks a bit mysterious: it's neither a straightforward variant of gradient descent, nor of mirror
descent, nor of proximal gradient descent.

It's strongly reminiscent, though, of the Frank-Wolfe a.k.a conditional
gradient method,[^3] since the primal update consists in setting
$$w_{t+1}$$ to a convex combination of $$w_t$$ and $$u_t$$. 
And indeed, it may be seen as a generalization of Frank-Wolfe to regularized instead of
constrained optimization problems ([Bach 2013](http://arxiv.org/abs/1211.6302), equation (17)).
To see this, consider the case where $$\Psi(w) = \iota_\Omega(w)$$ for some
convex set $$\Omega$$. Then, denoting $$g_t = V^* a_t$$, the above method is
equivalent to 

$$\begin{aligned}
    g_t &= V^* a_t 
    = V^* \nabla \mathcal{L}(V w_t) 
    = \left.\nabla_w \mathcal{L}(V w)\right|_{w_t} \\
    u_t &= \partial \Psi^*(-g_t)
    = \mathop{\mathrm{arg\,max}}_s \left\langle s, -g_t \right\rangle - \Psi(s)
    = \mathop{\mathrm{arg\,min}}_{s \in \Omega} \left\langle s, g_t \right\rangle \\
    w_{t+1} &= (1-\sigma_t) w_t + \sigma_t u_t\end{aligned}$$
    
which is exactly the Frank-Wolfe algorithm for the problem
$$\min_w \mathcal{L}(Vw) ~\text{s.t}~ w \in \Omega$$.

### Equivalence to a fully dual approach

Denote $$S_P$$ (resp. $$S_D$$) the optimal solution set for the primal
problem $$\eqref{eq:FRDT_primal}$$ (resp. dual problem
$$\eqref{eq:FRDT_dual}$$). According to the FRDT, $$w \in S_P$$ iff there
exists $$a \in S_D$$ such that
$$\eqref{eq:FRDT_KKT_L}$$. This motivates the following fully dual
approach:

-   Choose an update rule for $$w_{t+1}$$ such that
    $$\eqref{eq:FRDT_KKT_L}$$ is always satisfied;

-   Choose an update rule for $$a_{t+1}$$ that takes a step towards
    solving the dual problem: $$\max_{a'} D(a')$$ i.e
    $$\min_{a'} -D(a') = \Psi^*(-V^* a') + \mathcal{L}^*(a')$$.

For the dual update rule we may choose mirror descent using
$$\mathcal{L}^*$$.
For the primal update rule, we may apply the same trick as above to enforce
$$\eqref{eq:FRDT_KKT_L}$$
Thus we are led to the exact same algorithm as by the duality-gap approach.

In hindsight this equivalence is not surprising at all. Indeed, since we
enforce $$\eqref{eq:FRDT_KKT_L}$$, the dual problem (towards solving which we
let $$a_{t+1}$$ take a step) is the same in both approaches:

$$
\min_{a'} F(w_t,a')-D(a')
    = P(w_t)-D(a')
    \equiv \max_{a'} D(a').
$$

[^1]: <https://blogs.princeton.edu/imabandit/2013/04/16/orf523-mirror-descent-part-iii/>

[^2]: <https://blogs.princeton.edu/imabandit/2013/04/23/orf523-mirror-prox/>

[^3]: <https://en.wikipedia.org/wiki/Frank-Wolfe_algorithm>
