# Tensor
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
Class `measure` based on `object` class is providing unified logic for working with state variable.


# State
Class `state` has array of `measure` that describe a state of an object.
State has only one `basis` object so any tensor measures are linked to the same basis, at the same time indifferent scalar measures are not linked to any basis.
