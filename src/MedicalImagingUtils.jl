module MedicalImagingUtils

using Images

include("./calciumscoring.jl")

export
    # Export calciumscoring.jl functions
    compute_calcium_scores

end
