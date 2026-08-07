module BilinearProjections

using LinearAlgebra

export AbstractBilinearConstraint,
       CrossConstraint,
       BilinearLevelSet,
       CanonicalHyperbolicParaboloid,
       HyperbolicParaboloid,
       AbstractProjectionSet,
       SingletonProjection,
       SphereProjection,
       ProjectionInfo,
       AbstractSelector,
       DefaultSelector,
       UnitDirection,
       projection_set,
       project,
       project!,
       select,
       residual,
       isfeasible,
       distance2

include("types.jl")
include("utilities.jl")
include("root_solvers.jl")
include("canonical.jl")
include("interface.jl")
include("cross.jl")
include("level_set.jl")
include("paraboloid.jl")

end
