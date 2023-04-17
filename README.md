# Tensor
Class `object` is a handler of invariant tensor object that is represented by components and basis at any time.

Base class `container` may contain of arrays with any dimension and rank (even rank = 0 for scalar arrays).
Based on `container` class `object` provides basic tensor-vector and scalar arrays calculus in N dimensional spatial.

Class `object` provide +/-/scal product operators between each other. So you should not keep in mind at how basis are component of current tensor/vector, 
internal functionality does operation correctly.

Using smart pointers for manage basis and component containers allows to have many objects (tensor/vector) at the same basis object. 
So all objects linked to the one basis will change at the same time when basis will be changed.
Basis object will live till the last tensor/vector object will be destroyed.

There is possibility to handle array of scalar or even just scalar variable like tensor object wut without basis poiter.
It is an advantage for abstraction and construction of any measures (see the next paragraph).

# Measure
The main purposes the next features are imtroduced for fast implementation of any mathematical model based on [state variable approach](https://en.wikipedia.org/wiki/State_variable). The main points:
- state of a system may be described by finite set of tensor or/and scalar state variables (**SV**)
- each state variable may have evolutionary equation (usually - differential equation)
- each state variable may have dependecies on other strate variables
- 
Class `StateMeasure` based on `object` class is providing unified logic for working with state variable:
- it has `rate` - value of rate of **SV**
- it has `value` - value **SV**
- also it contains previous value and rate of **SV** for ability to implement of more convinience numerical schemas 
- `rate_equation()` - evolution equation in rate form (the rate of state variable dependency on set of parameters and set of **SV**'s)
- `finit_equation()` - evolution equation in finite form (describes dependency of finite value on set of parameters and set of **SV**'s)
- `calc_rate()` - by default it is the first order schema to calculate rate, may be overriden
- `integrate_value(T dt)` - by default the first order (Euler) schema to integrate value, may be overriden

# Implemented measures
There are a few predefined wide used strain and stress measures:
- deformation gradient `GradDeform`, it is a base measure in a [Finite strain therory](https://en.wikipedia.org/wiki/Finite_strain_theory)
- Cauchy stress measure `CaushyStress`, it is a base measure [Caushy stress](https://en.wikipedia.org/wiki/Cauchy_stress_tensor) from wich (and `GradDeform`) my be derivated another stress measures like Piola–Kirchhoff tensor.

# Numerical schema
Base abstract class `AbstractSchema` describes main functionality:
- `init()` - assumed that it is called before calc rate/finite and integration procedure
- `calc()` - the main stage of numerical schema, `rate` or `value` on the end of current step are calculated
- `finalize()` - actions after if necessery

There are a few types of numerical schema:
- `RATE_CALCULATE`: *dX(n+1) := F(...), X(n+1) = X(n) + dX(n+1)dt*
- `FINITE_CALCULATE`: *X(n+1) := G(...), dX(n+1) = (X(n+1)-X(n))/dt*

Class `DefaultSchema` implements (abstract methods of `AbstractSchema`) a plain first order numerical schema over `StateMeasure`.

# State
Class `State` has array of `measure`'s that describes a state of an object.
State has only one `basis` object so any tensor measures are linked to the same basis, at the same time indifferent scalar measures are not linked to any basis.

# Models
Any class inherited from class `State` is a material model contains array of `measure`'s, logic of any neccesary calculation and `Relation` links two `measure`'s. 

# Relation
Class `Relation` implements relationship between two arbitrary object inherited from `DefaultSchema`. All relation have rate and finite implementation and may be used in mix.

## Base Elastic Relation
`ElasticRelation` inherited from `Relation` and implements relationship between `CaushyStress` and `GradDeform` in any combination.
There are two relations:
- `HypoElasticRelation` inherited from `ElasticRelation` and implements relationship between `CaushyStress` and lagrangian strain tensor (derived from `GradDeform`)
- `HyperElasticRelation` inherited from `ElasticRelation` and implements relationship between `CaushyStress` and spatial velocity gradient (derived from rate of `GradDeform`)

## Inelastic (plastic) Relation
In fact there is no specific Relation like `HypoElasticRelation`/`HyperElasticRelation` due to describe plastic behavior. Instead of it there is used `GradDeformInelast` measure that describes how Plastic part will be calculated.

## Strain decomposition
The same as for Inelastic Relation. There is used `GradDeformElast`measure that describes how elastic part is extracted from full and plastic measure.
