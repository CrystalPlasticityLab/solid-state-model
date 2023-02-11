# tensor
Classes Tensor and vector provides basic tensor-vector calculus in N dimensional spatial

Classes Tensor and vector provide +/-/scal product operators between each other. So you should not keep in mind at how basis are component of current tensor/vector, 
internal functionality does operation correctly.

Using shared_ptr for manage basis object allows to have many objects (Tensor/vector) at the same basis object. 
So all objects linked to the one basis will change at the same time when basis will be changed.
Basis object will live till the last Tensor/vector object will be destroyed.
