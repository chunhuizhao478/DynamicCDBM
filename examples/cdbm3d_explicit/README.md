# 3D Continuum Damage-Breakage Model Explicit Time Integration Prototype
## Created By Chunhui Zhao, Aug 11th

### Highlights:
- The fault is modeled as high damage strip 
- The fault is buried in the medium
- The nucleation is done by set a over-damaged patch inside the strip
- The initial condition is obtained by (AD) static solve with boundary traction applied to the domain
- Explicit Time Integration (CentralDifference) is used for dynamic solve
- NonReflecting Boundary Condition is applied on all boundary surfaces