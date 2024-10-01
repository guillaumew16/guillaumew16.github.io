---
layout: post
date:   2024-10-01
categories: math
author: Guillaume Wang
title: Two convergence analyses of "self"-mirror flow
---

The mirror descent algorithm, in its simplest form, is given by the
update rule 

$$\nabla h(x^{k+1}) = \nabla h(x^k) - \eta \nabla f(x^k),$$

where $$f$$ is the objective function, $$\eta>0$$ is the step-size, and $$h$$
is a hyperparameter of the algorithm called the link function, such that
$$\nabla h: \mathbb{R}^d \to \mathbb{R}^d$$ is bijective. The convergence
analysis of mirror descent (i.e., the rate at which $$f(x^k) \to \inf f$$
as $$k \to \infty$$) is a well-studied question, dating back at least to
Nemirovsky and Yudin who first introduced it [1].
Sebastian Bubeck's book has a nice summary of the results
[2, Chapter 4]. In a more recent development (2017), the
works [3,4] offered a fresh
perspective on the question, which is the one I will use.

In the limit where $$\eta \to 0$$, the sequence $$(x^k)_{k \geq 0}$$ defines
a continuous curve via $$x(t) = \lim_\eta x^{\lfloor{t/\eta}\rfloor}$$, called the
mirror flow. (This is a generalization of gradient flow, which is the
continuous limit of gradient descent.) AFAIK, its first explicit
appearance in the literature is the work by Amid and Warmuth
[7]. The convergence analysis of mirror flow has
not been written down explicitly in a research paper, but it follows
easily from adapting the proofs of [4].

The aforementioned stuff comes from the optimization literature. More or
less independently, people interested in sampling and/or PDEs have been
looking at the "birth-death dynamics" equation
[5,6,9]

$$\forall i \in \{1,...,N\},~ \frac{d}{dt} \mu_t[i] = -\mu_t[i] \left( \log \frac{\mu_t[i]}{\nu[i]} - \mathsf{KL}\left( \mu_t \middle\| \nu \right) \right),$$

where
$$\mu_t \in \Delta_N = \left\{ \mu \in \mathbb{R}_+^N, \sum_i \mu[i] = 1 \right\}$$
(the probability simplex), $$\nu$$ is a fixed element of $$\Delta_N$$, and
$$\mathsf{KL}\left( \mu \middle\| \nu \right) = \sum_j \mu[j] \log \frac{\mu_t[i]}{\nu[i]}$$.
[^1] The KL term is there to ensure $$\mu_t$$ remains normalized. As we
will see, this equation can be interpreted as constrained mirror flow of
$$f(\mu) = \mathsf{KL}\left( \mu \middle\| \nu \right)$$, with link
function
$$h(\mu) = \sum_i \mu[i] \log \mu[i] = \mathsf{KL}\left( \mu \middle\| {\boldsymbol{1}} \right)$$.

If we apply the results from the optimization literature, we get the
convergence guarantees

$$\mathsf{KL}\left( \nu \middle\| \mu_t \right) \leq e^{-t} \mathsf{KL}\left( \nu \middle\| \mu_0 \right)
~~~~\text{and}~~~~
\mathsf{KL}\left( \mu_t \middle\| \nu \right) \leq \frac{1}{e^t - 1} \mathsf{KL}\left( \nu \middle\| \mu_0 \right),$$

see [10] for example. I was surprised to learn that the
sampling/PDE community obtained convergence guarantees of a very
different form: in [9],
$$\mathsf{KL}\left( \mu_t \middle\| \nu \right) \leq e^{-t} \cdot \mathrm{cst}( \min_i \frac{\mu_0[i]}{\nu[i]}, \mathsf{KL}\left( \mu_0 \middle\| \nu \right))$$,
and in [8], the sharper result
$$\mathsf{KL}\left( \mu_t \middle\| \nu \right) = \frac{\kappa_2}{2} e^{-2t} + O(e^{-3t})$$
for some constant $$\kappa_2$$---yes, this is an equality!---and even more
precisely, 

$$\label{eq:KL_expan} \tag{1}
\mathsf{KL}\left( \mu_t \middle\| \nu \right) = \frac{\kappa_2}{2} e^{-2t} + \sum_{n=3}^\infty \frac{\kappa_n}{n(n-2)!} e^{-nt}$$

where
$$\kappa_n = {\left.\frac{d^n}{dz^n}\right|_{z=0}} \log \sum_i \nu[i] \exp\left( z \log \frac{\mu_0[i]}{\nu[i]} \right)$$.

In this blog post,

-   I review the convergence analysis of mirror flow following
[4]. The only quantitative assumption is that $$f$$ is
relatively strongly convex w.r.t. $$h$$.

-   I show that an expansion analogous to
$$\eqref{eq:KL_expan}$$ holds for any "self"-mirror flow, i.e.,
when $$f = h + \mathrm{linear}$$.

Of course, similar things hold for the mirror descent algorithm, rather
than mirror flow.

Thanks go to Carles and Aram (the authors of [8]) for
interesting discussions. They were aware of the fact that the expansion
$$\eqref{eq:KL_expan}$$ could be generalized to all self-mirror flows,
but decided it was of limited research value; and I agree with them, but
it's the kind of fun fact that's worth sharing.

# Mirror flow (with linear constraint set) definitions

In this post, we are interested in a (constrained) optimization problem
of the form 

$$\label{eq:pb} \tag{2}
\min_{x \in \mathcal{Z}} f(x)$$

where the constraint set
$$\mathcal{Z}$$ is linear, i.e.,
$$\mathcal{Z}= \left\{ x \in \mathbb{R}^d ; Ax = b \right\}$$ for some
$$A \in \mathbb{R}^{m \times d}$$ and $$b \in \mathbb{R}^m$$. Furthermore,
we assume that we have first-order access to $$f$$, as well as to $$h$$ and
$$h^*$$ for some strictly convex "link function"
$$h: \mathcal{Z}\to \mathbb{R}$$. For simplicity we assume $$f$$, $$h$$ and
$$h^*$$ are $$C^3$$.

<div class="definition" markdown="1" text="[7]">
For an initial point
$$x_0 \in \mathcal{Z}$$, the Mirror Flow (MF) of $$f$$ with link function
$$h$$ is the unique curve $$(x_t)_t$$ such that $$x_{t=0} = x_0$$ and

$$\frac{dx_t}{dt} = -\Phi_{x_t}^{-1} P_{x_t} \nabla f(x_t)$$ 

where
$$\Phi_x = \nabla^2 h(x)$$ and
$$P_x = I_d - A^\top \left[ A \Phi_x^{-1} A^\top \right]^{-1} A \Phi_x^{-1}$$.
</div>

<div class="remark" markdown="1">
*Remark 1*. One can check that $$P_x^2 = P_x$$ and
$$P_x \Phi_x = \Phi_x P_x^\top$$, and that $$P_x^\top$$ is the orthogonal
projection onto $$\mathop{\mathrm{Ker}}A$$ w.r.t. the inner product
$$\left\langle \cdot, \cdot \right\rangle_{\Phi_x}$$. [^2]
</div>

When analyzing the convergence of MF, the following conditions turn out
to be quite natural.

<div class="definition" markdown="1" text="[4]">
Let
$$\mathcal{Z}\subset \mathbb{R}^d$$ be a convex set,
$$f, h: \mathcal{Z}\to \mathbb{R}\cup \{\infty\}$$ and $$L, \mu \geq 0$$.
Denote by $$\mathrm{relint}(\mathcal{Z})$$ the relative interior of
$$\mathcal{Z}$$ and by
$$D_f(x, x') = f(x) - f(x') - \nabla f(x')^\top (x-x')$$ the Bregman
divergence of $$f$$.

$$f$$ is called $$L$$-smooth relatively to $$h$$ if
$$\forall x, x' \in \mathrm{relint}(\mathcal{Z}), D_f(x, x') \leq L D_h(x, x')$$.

$$f$$ is called $$\mu$$-strongly convex relatively to $$h$$ if
$$\forall x, x' \in \mathrm{relint}(\mathcal{Z}), D_f(x, x') \geq \mu D_h(x, x')$$.
</div>

In this post, we are interested in the convergence of MF
for $$\eqref{eq:pb}$$ under
the assumption of relative strong convexity:

<div class="assumption" markdown="1">
$$f$$ is $$\mu$$-strongly convex relatively to $$h$$.
Moreover, $$f$$ has a unique minimizer, and we denote
$$x^* = \mathop{\mathrm{arg\,min}}f$$ and $$f^* = f(x^*)$$. We also denote
by $$(x_t)_t$$ a MF of $$f$$ with link function $$h$$.
</div>

<div class="remark" markdown="1">
We have $$f(x) - f^* = D_f(x, x^*)$$. Indeed,

$$f(x) - f^* - D_f(x, x^*) 
= \nabla f(x^*)^\top (x - x^*) 
= \Big( P_{x^*} \nabla f(x^*) \Big)^\top (x - x^*) 
= 0.$$ 

The second equality uses that
$$P_{x^*}^\top (x - x^*) = x-x^*$$ since
$$x-x^* \in \mathop{\mathrm{Ker}}A$$, and the third equality uses that
$$\Phi_{x^*}^{-1} P_{x^*} \nabla f(x^*) = 0$$ since $$x^*$$ is a stationary
point of MF.
</div>

# The classical optimization analysis

The following proposition establishes contraction in $$D_h(x^*, \cdot)$$.
It appeared as [10, Theorem 1] in the context of
Fisher-Rao gradient flows (generalizations of birth-death dynamics).

<div class="proposition" markdown="1">
We have
$$D_h(x^*, x_t) \leq e^{-\mu t} D_h(x^*, x_0)$$.
</div>

<div class="proof" markdown="1">
We have 

$$\frac{d}{dt} D_h(x^*, x_t) 
= \nabla f(x_t)^\top (x^*-x_t)
= f^* - f(x_t) - D_f(x^*, x_t)
\leq -D_f(x^*, x_t)
\leq -\mu D_h(x^*, x_t),$$ 

and the proposition follows by
applying Grönwall's lemma. The only non-obvious step in this string of
(in)equalities is the first one, i.e., the fact that
$$\frac{d}{dt} D_h(x^*, x_t) = \nabla f(x_t)^\top (x^*-x_t)$$. For this,
just substitute the definition and compute:

$$\begin{aligned}
\frac{d}{dt} D_h(x^*, x_t) 
&= \frac{d}{dt} \left[ h(x^*) - h(x_t) - \nabla h(x_t)^\top (x^*-x_t) \right] \\
&= 0 - \nabla h(x_t)^\top \dot{x}_t - \nabla h(x_t)^\top (0-\dot{x}_t)
- \dot{x}_t^\top \nabla^2 h(x_t) (x^* - x_t) \\
&= 0 + \Big( \Phi_{x_t}^{-1} P_{x_t} \nabla f(x_t) \Big)^\top \nabla^2 h(x_t) (x^* - x_t) \\
&= \nabla f(x_t)^\top P_{x_t}^\top (x^* - x_t)
~~~~ ~~\text{since}~~ \Phi_x = \nabla^2 h(x).
\end{aligned}$$ 

To conclude, note that
$$P_{x_t}^\top (x^* - x_t) = x^* - x_t$$ since
$$x^* - x_t \in \mathop{\mathrm{Ker}}A$$ since
$$x^*, x_t \in \mathcal{Z}$$.
</div>

The above result can be deduced by taking the small-timestep limit of
the original argument of [4, Theorem 3.1]; in fact by
adapting more attentively the original argument, we have the following
convergence in function value.

<div class="proposition" markdown="1">
We have that $$t \mapsto f(x_t)$$ is non-increasing and 

$$f(x_t) - f^* \leq \frac{\mu}{e^{\mu t}-1}
% \left( D_h(x^*, x_0) - e^{\mu t} D_h(x^*, x_t) \right).
D_h(x^*, x_0).$$
</div>


<div class="remark" markdown="1">
Observe that
$$\lim_{\mu \to 0} \frac{\mu}{e^{\mu t}-1} = \big( {\left.\frac{d}{d\mu}\right|_{\mu=0}} e^{\mu t} \big)^{-1} = \frac1t$$.
This suggests that when $$f$$ is only convex, we still have convergence of
MF in function value with the rate $$1/t$$. This is indeed the case, as
could be shown by adapting the arguments.
</div>

<div class="proof" markdown="1">
By definition of MF,

$$\frac{d}{dt} f(x_t) = -\nabla f(x_t)^\top \Phi_{x_t}^{-1} P_{x_t} \nabla f(x_t).$$

Since $$P_x^2 = P_x$$ and $$P_x \Phi_x = \Phi_x P_x^\top$$, then
$$\Phi_x^{-1} P_x = (\Phi_x^{-1} P_x) P_x = P_x^\top \Phi_x^{-1} P_x$$ is
symmetric positive-semidefinite. Hence $$\frac{d}{dt} f(x_t) \leq 0$$.

The previous proposition showed that 

$$\begin{aligned}
\forall s,~ \frac{d}{ds} D_h(x^*, x_s)
&= f^* - f(x_s) - D_f(x^*, x_s) \\
f(x_s) - f^* &= -D_f(x^*, x_s) - \frac{d}{ds} D_h(x^*, x_s) \\
&\leq -\mu D_h(x^*, x_s) - \frac{d}{ds} D_h(x^*, x_s) \\
&= -e^{-\mu s} \frac{d}{ds} \left[ e^{\mu s} D_h(x^*, x_s) \right] \\
\text{and so}~~
\int_0^t e^{\mu s} \left( f(x_s) - f^* \right) ds
&\leq -\int_0^t \frac{d}{ds} \left[ e^{\mu s} D_h(x^*, x_s) \right] ds \\
&= D_h(x^*, x_0) - e^{\mu t} D_h(x^*, x_t). 
\end{aligned}$$

On the other hand, since $$t \mapsto f(x_t)$$ is non-increasing,

$$\int_0^t e^{\mu s} \left( f(x_s) - f^* \right) ds
\geq \left( f(x_t) - f^* \right) \int_0^t e^{\mu s} ds
= \left( f(x_t) - f^* \right) \frac1\mu (e^{\mu t}-1).$$

The proposition follows by combining the two inequalities above. 
</div>

# Self-mirror flow

Let us now consider the convergence behavior of MF for
$$\eqref{eq:pb}$$ when $$f$$
is equal to $$h$$ itself up to affine terms. In particular,
$$D_f(x, x') = D_h(x, x')$$, and so $$f$$ is of course $$1$$-strongly convex
relatively to $$h$$. In fact we can write

$$f(x) = D_f(x, x^*) + f^* = D_h(x, x^*) + f^*$$ 

so we will assume w.l.o.g. that
$$f = D_h(\cdot, x^*)$$.

In this case, the convergence guarantees from the classical analysis
become 

$$\frac{d}{dt} D_f(x^*, x_t)
= f^* - f(x_t) - D_f(x^*, x_t)
\leq -D_f(x^*, x_t)$$ 

and in particular
$$D_f(x^*, x_t) \leq e^{-t} D_f(x^*, x_0)$$, and 

$$\begin{aligned}
f(x_t) - f^* 
\leq \frac{1}{e^t-1} \int_0^t e^s \left( f(x_s) - f^* \right) ds
= \frac{1}{e^t-1} \left( D_f(x^*, x_0) - e^t D_f(x^*, x_t) \right).\end{aligned}$$

But that still only gives us upper bounds. We actually have much better:
the MF trajectory $$(x_t)_t$$ coincides, up to time-reparametrization,
with straight lines in the dual space. This is quite intuitive in the
unconstrained case, and for linear equality constraints the only
potential difficulty is to properly specify what is meant by "dual".

#### Unconstrained case.

One can check by applying the definitions the following.

<div class="proposition" markdown="1">
Suppose $$\mathcal{Z}= \mathbb{R}^d$$. Then MF is
equivalent to 

$$x_t = \nabla h^*(y_t)
~~\text{and}~~
\frac{dy_t}{dt} = -\nabla f(x_t)$$ 

where $$h^*$$ denotes the
convex conjugate of $$h$$, characterized by
$$\nabla h^* = (\nabla h)^{-1}$$.

Suppose furthermore that $$f(x) = D_h(x, x^*)$$ for some $$x^*$$, and let
$$y^* = \nabla h(x^*)$$. Then $$\nabla f = \nabla h - y^*$$ and
$$\frac{dy_t}{dt} = y^* - \nabla h(x_t) = -(y_t - y^*)$$, so
$$y_t = y^* + e^{-t} (y_0 - y^*)$$.
</div>

#### With linear equality constraints.

Recall that the link function $$h$$ is defined as a function
$$\mathcal{Z}\to \mathbb{R}$$, so the notation $$\nabla h^*$$ itself
deserves some clarification.

<div class="definition" markdown="1">
The gradient of the strictly convex and differentiable
function $$h: \mathcal{Z}\to \mathbb{R}$$ is the mapping
$$\nabla h: \mathcal{Z}\to \mathop{\mathrm{Ker}}A$$ such that

$$\forall v \in \mathop{\mathrm{Ker}}A,~
h(x + \varepsilon v) - h(x) = \varepsilon v^\top \nabla h(x) + o(\varepsilon).$$

The convex conjugate of $$h$$ is the function
$$h^*: \mathop{\mathrm{Ker}}A \to \mathbb{R}$$ defined by
$$h^*(y) = \sup_{x \in \mathcal{Z}} x^\top y - h(x)$$, and $$\nabla h^*$$ is
defined similarly to $$\nabla h$$.

One can show that it holds
$$(\nabla h^*) \circ (\nabla h) = \mathop{\mathrm{id}}_{\mathcal{Z}}$$ and
$$(\nabla h) \circ (\nabla h^*) = \mathop{\mathrm{id}}_{\mathop{\mathrm{Ker}}A}$$.
</div>

<div class="proposition" markdown="1">
Denote by $$\Pi$$ the orthogonal projector from
$$\mathbb{R}^d$$ onto $$\mathop{\mathrm{Ker}}A$$. MF is equivalent to
$$x_t = \nabla h^*(y_t)$$, i.e.,
$$y_t = \nabla h(x_t) \in \mathop{\mathrm{Ker}}A$$, and 

$$\frac{dy_t}{dt} 
= -\nabla^2 h(x_t) \Phi_{x_t}^{-1} P_{x_t} \nabla f(x_t)
= -P_{x_t} \nabla f(x_t)
= \Pi \frac{dy_t}{dt}
= -\Pi P_{x_t} \nabla f(x_t)
= -\Pi \nabla f(x_t).$$

Suppose furthermore that $$f(x) = D_h(x, x^*)$$ for some $$x^*$$ and let
$$y^* := \nabla h(x^*) \in \mathop{\mathrm{Ker}}A$$. Then
$$\nabla f = \nabla h - y^*$$ and
$$\frac{dy_t}{dt} = -\Pi\left( \nabla h(x_t) - y^* \right) = -(y_t - y^*)$$,
so $$y_t = y^* + e^{-t} (y_0 - y^*)$$.
</div>

<div class="proof" markdown="1">
The first part of the proposition follows from the definition
of MF and from the fact that
$$\Pi P_x = \Pi (I_d - A^\top \left[ A \Phi_x^{-1} A^\top \right]^{-1} A \Phi_x^{-1}) = \Pi$$,
since $$A \Pi^\top = A \Pi = 0$$. The rest of the proposition can be
checked directly. ◻
</div>

#### The expansion.

Consider the unconstrained case, and suppose $$f(x) = D_h(x, x^*)$$. In
particular $$\nabla f = \nabla h - y^*$$, or equivalently,
$$\nabla f^* = \nabla h^*(\cdot + y^*)$$. We have the closed-form formula,
for all $$t \geq 0$$, 

$$x_t = \nabla h^*(y_t)
= \nabla h^*(y^* + e^{-t} (y_0 - y^*))
= \nabla f^*(e^{-t} (y_0 - y^*)).$$ 

So let us pose $$Y = y_0 - y^*$$
and, for $$\tau \in [0,1]$$, $$Y_\tau = (1-\tau) Y$$ and
$$\tilde{x}_\tau = \nabla f^*(Y_\tau)$$. Then $$\tilde{x}_0 = x_0$$ and
$$\tilde{x}_1 = x^*$$, and by the equality case in Fenchel-Young's
inequality, 

$$\begin{aligned}
f(\tilde{x}_\tau)
= f(\nabla f^*(Y_\tau))
&= -f^*(Y_\tau) + Y_\tau^\top \nabla f^*(Y_\tau) \\
&= -\sum_{n=0}^\infty \frac{1}{n!} \nabla^n f^*(0) : Y_\tau^{\otimes n}
+ Y_\tau^\top \left[ \sum_{n=0}^\infty \frac{1}{n!} \nabla^{n+1} f^*(0) \cdot Y_\tau^{\otimes n} \right]\\
&= -f^*(0) + \sum_{n=1}^\infty \left[ -\frac{1}{n!} + \frac{1}{(n-1)!} \right] \nabla^n f^*(0) : Y_\tau^{\otimes n} \\
&= (\inf f) + \sum_{n=2}^\infty \frac{1}{n (n-2)!} \left[ \nabla^n f^*(0) : Y^{\otimes n} \right] (1-\tau)^n.\end{aligned}$$

Hence, substituting $$t$$ such that $$1-\tau = e^{-t}$$, we have
$$\tilde{x}_\tau = x_t$$ and

$$f(x_t) - (\inf f) = \sum_{n=2}^\infty \frac{1}{n (n-2)!} \left[ \nabla^n f^*(0) : Y^{\otimes n} \right] e^{-nt}.$$

In the case with linear equality constraints, actually the same
derivation applies verbatim. So we have proved:

<div class="theorem" markdown="1">
Suppose $$f(x) = D_h(x,x^*)$$ for some
$$x^* \in \mathcal{Z}$$. Denote by $$f^*$$ the convex conjugate
of $${\left.f\right|_{\mathcal{Z}}}$$, defined by
$$f^*(y) = \sup_{x \in \mathcal{Z}} x^\top y - f(x)$$. Then

$$f(x_t) = \sum_{n=2}^\infty \frac{1}{n (n-2)!} \left[ {\left.\frac{d^n}{dz^n}\right|_{z=0}} f^*(zY) \right] e^{-nt} 
~~~~\text{where}~~
Y = \nabla h(x_0) - \nabla h(x^*).$$
</div>

#### Birth-death dynamics.

To recover the result of [8], consider
$$\mathcal{Z}= \Delta_N = \left\{ \mu \in \mathbb{R}_+^N; {\boldsymbol{1}}^\top \mu = 1 \right\}$$
and $$f, h: \mathcal{Z}\to \mathbb{R}$$ defined by
$$f(\mu) = \mathsf{KL}\left( \mu \middle\| \nu \right)$$ and
$$h(\mu) = \mathsf{KL}\left( \mu \middle\| {\boldsymbol{1}} \right)$$. 
One can check that MF for this $$f$$ with this link function $$h$$, is precisely the birth-death dynamics from the introduction.

It is well-known that $$f(\mu) = D_h(\mu, \nu)$$ and that

$$f^*(y) = \sup_{\mu \in \Delta_N} \sum_i \mu[i]~ y[i] - \mathsf{KL}\left( \mu \middle\| \nu \right)
= \log \sum_i \nu[i]~ \exp(y[i])$$ 

(this identity is sometimes called Donsker-Varadhan duality). Moreover, $$Y[i] := \big( \nabla h(\mu_0) - \nabla h(\nu) \big)[i] = \log \frac{\mu_0[i]}{\nu[i]}$$. So the theorem applied to this $$f$$ and $$h$$ asserts that, for $$(\mu_t)_t$$ following the birth-death dynamics,

$$f(\mu_t) = \sum_{n=2}^\infty \frac{1}{n (n-2)!} ~\underbrace{
\left[ {\left.\frac{d^n}{dz^n}\right|_{z=0}} \log \sum_i \nu[i]~ \exp\left( z  \log \frac{\mu_0[i]}{\nu[i]} \right) \right]
}_{\kappa_n}~
e^{-nt},$$

which is precisely the statement of the main result of [8].

---


**References**

[1] Arkadij Semenovič Nemirovskij and David Borisovich Yudin. “Problem complexity and method efficiency in optimization”. In: (1983).

[2] Sébastien Bubeck. “Convex optimization: Algorithms and complexity”. In: Foundations and Trends in Machine Learning 8.3-4 (2015), pp. 231–357.

[3] Heinz H Bauschke, Jérôme Bolte, and Marc Teboulle. “A descent lemma beyond Lipschitz gradient continuity: first-order methods revisited and applications”. In: Mathematics of Operations Research 42.2 (2017), pp. 330–348.

[4] Haihao Lu, Robert M Freund, and Yurii Nesterov. “Relatively smooth convex optimization by first-order methods, and applications”. In: SIAM Journal on Optimization 28.1 (2018), pp. 333–354.

[5] Yulong Lu, Jianfeng Lu, and James Nolen. “Accelerating langevin sampling with birth-death”. In: arXiv preprint arXiv:1905.09863 (2019).

[6] Grant Rotskoff, Samy Jelassi, Joan Bruna, and Eric Vanden-Eijnden. “Global convergence of neuron birth-death dynamics”. In: International Conference on Machine Learning. 2019.

[7] Ehsan Amid and Manfred KK Warmuth. “Reparameterizing mirror descent as gradient descent”. In: Advances in Neural Information Processing Systems 33 (2020), pp. 8430–8439.

[8] Carles Domingo-Enrich and Aram-Alexandre Pooladian. “An explicit expansion of the Kullback-Leibler divergence along its Fisher-Rao gradient flow”. In: Transactions on Machine Learning Research (2023).

[9] Yulong Lu, Dejan Slepčev, and Lihan Wang. “Birth–death dynamics for sampling: global convergence, approximations and their asymptotics”. In: Nonlinearity 36.11 (2023), p. 5731.

[10] Rentian Yao, Linjun Huang, and Yun Yang. “Minimizing Convex Functionals over Space of Probability Measures via KL Divergence Gradient Flow”. In: International Conference on Artificial Intelligence and Statistics. PMLR. 2024, pp. 2530–2538.

------

[^1]: Actually people often look at a "continuous-space" version of this equation, where $$(\mu_t[i])_i \in \Delta_N$$ is replaced by $$(\mu_t(x))_{x \in \mathbb{R}^d}$$ in the space of probability *measures*, but the convergence behavior is the same except for regularity issues which can be treated separately.

[^2]: To check this, note that $$P_x A^\top = 0$$ and $$\mathop{\mathrm{Ker}}P_x = \mathop{\mathrm{Im}}(I_d - P_x) \subset \mathop{\mathrm{Im}}A^\top$$.
