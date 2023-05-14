# Tensor
Class `object` is a handler of invariant tensor object that is represented by components and basis at any time.

Base class `container` may contain of arrays with any dimension and rank (even rank = 0 for scalar arrays). 
Class `container` is inherited from `std::array`, so it is fast.
Based on `container` class `object` provides basic tensor-vector and scalar arrays calculus in N dimensional spatial.

Class `object` provide +/-/scal product operators between each other. So you should not keep in mind at how basis are component of current tensor/vector, 
internal functionality does operation correctly.

Using smart pointers for manage basis and component containers allows to have many objects (tensor/vector) at the same basis object. 
So all objects linked to the one basis will change at the same time when basis will be changed.
Basis object will live till the last tensor/vector object will be destroyed.

There is possibility to handle array of scalar or even just scalar variable like tensor object but without basis poiter.
It is an advantage for abstraction and construction of any measures (see the next paragraph).

# Measure
The main purposes the next features are imtroduced for fast implementation of any mathematical model based on [state variable approach](https://en.wikipedia.org/wiki/State_variable). The main points:
- state of a system may be described by finite set of tensor or/and scalar state variables (**SV**)
- each state variable may have evolutionary equation (usually - differential equation)
- each state variable may have dependecies on other strate variables

Class `StateMeasure` based on `object` class is providing unified logic for working with state variable:
- it has `rate` - value of rate of **SV**
- it has `value` - value **SV**
- also it contains previous value and rate of **SV** for ability to implement of more convinience numerical schemas 
- `rate_equation(T t, T dt)` - evolution equation in rate form (the rate of state variable dependency on set of parameters and set of **SV**'s)
- `finit_equation(T t, T dt)` - evolution equation in finite form (describes dependency of finite value on set of parameters and set of **SV**'s)
- `calc_rate(T dt)` - by default it is the first order schema to calculate rate, may be overriden
- `integrate_value(T dt)` - by default the first order (Euler) schema to integrate value, may be overriden

# Implemented measures
There are a few predefined wide used strain and stress measures:
- deformation gradient `GradDeform`, it is a base measure in a [Finite strain therory](https://en.wikipedia.org/wiki/Finite_strain_theory), it also contains methods to calculate derived measures (left, right Hencky measures, stretch and Cauchy tensors and etc)
- Cauchy stress measure `CaushyStress`, it is a base measure [Caushy stress](https://en.wikipedia.org/wiki/Cauchy_stress_tensor) from wich (and `GradDeform`) my be derivated another stress measures like Piola–Kirchhoff tensor.

# Numerical schema
Base abstract class `AbstractSchema` describes main functionality:
- `init()` - assumed that it is called before calc rate/finite and integration procedure
- `calc()` - the main stage of numerical schema, `rate` or `value` on the end of current step are calculated
- `finalize()` - actions after if necessery

There are a few types of numerical schema:
- `RATE_CALCULATE`: *dX(n+1) := F(...), X(n+1) = X(n) + dX(n+1)dt*
- `FINITE_CALCULATE`: *X(n+1) := G(...), dX(n+1) = (X(n+1)-X(n))/dt*

Class `StateMeasureSchema` implements (abstract methods of `AbstractSchema`) a plain first order numerical schema over `StateMeasure`.

# Material point

Class `MaterialPoint` is base for any material model. it is inherited from `AbstractSchema` and contains `basis` object (common for any measures), so any tensor measures are linked to the same basis, at the same time indifferent scalar measures are not linked to any basis.

# Models
Any class inherited from class `MaterialPoint` is a material model contains array of `StateMeasureSchema`'s, logic of any neccesary calculation and `Relation` links two or more `StateMeasureSchema`'s. 

# Relation
In case of dependency on `StateMeasureSchema` from other (or others) derived classes are implemeted:
- `ElasticRelation` (inherited from `StressMeasure`, template from `StressMeasure` and `StrainMeasure`) implemets elastic relation between `StressMeasure` and `StrainMeasure` (for example Hooke's law). It is extended stress measure of the material model (so it has polymorphic behavior as stress measure and relation).
- `PlasticRelation` (inherited from `StrainMeasure`, template parameters - `StressMeasure` and `StrainMeasure`) implements plastic flow rule (dependecy plastic strain part on stress). It is extended plastic strain measure of the material model (so it has polymorphic behavior as plastic strain measure and plastic relation).
- `StrainDecomposition` (inherited from `StrainMeasure`, template parameter - `StrainMeasure`) implements full strain measure decomposition into plastic and elastic ones. 
- `GeometricNonlinearity` implements dependency the basis orientation of `MaterialPoint` on the state. It may be a part of any model if necesessary.

## Base Elastic Model
`Elasticity` inherited from `MaterialPoint` (template parameters - `StressMeasure` and `StrainMeasure`) and implements elastic behavior. It has extended `ElasticRelation` instead of pure `StressMeasure` (it is both stress measure and relation with dependency on strain measure).

## Inelastic (plastic) Model
Class `Plasticity` inherited from `Elasticity` (template parameters - `StressMeasure` and `StrainMeasure`) and implements plastic behavior. It has `PlasticRelation` instead of pure plastic part of `Strain<easure` (it is both strain measure and relation with dependency on stress measure). Also it has `StrainDecomposition` as a rule for an elastic part calculation (full and plastic parts are known).

