### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 47d5b09c-c213-11eb-2472-7b23859c3caa
begin
	let
		using Pkg
		Pkg.activate(mktempdir())
		Pkg.Registry.update()
		Pkg.add("DICOM")
		Pkg.add("NIfTI")
		Pkg.add("Plots")
		Pkg.add(Pkg.PackageSpec(url="https://github.com/Dale-Black/DICOMUtils.jl"))
		Pkg.add("Images")
	end

	using PlutoUI
	using DICOM
	using NIfTI
	using Plots
	using DICOMUtils
	using LinearAlgebra
	using Images
end

# ╔═╡ e991c00b-974b-442c-93c6-55a0b12b27ed
TableOfContents()

# ╔═╡ 718b9504-f15f-41d1-b0ea-3346f17b4cca
md"""
## Load DICOMs
"""

# ╔═╡ a45f05b9-08d6-4277-9c9c-205f60b43391
volpath = "/Users/daleblack/Desktop/80.0";

# ╔═╡ 5e69faec-3f89-4f65-aab1-9429ad1e141b
volseg1 = "/Users/daleblack/Desktop/HEL_SLICER_SEG_0/S_1.20.nii";

# ╔═╡ 0399cca4-b251-451e-bf8c-96ccd207a1df
volseg2 = "/Users/daleblack/Desktop/HEL_SLICER_SEG_0/M_3.0.nii";

# ╔═╡ 20b7a9f9-4308-427f-995f-1f29f7f6df23
md"""
## Check Orientation
"""

# ╔═╡ e36c897d-07ae-4f61-8a5d-6dbdc651431a
vol = dcmdir_parse(volpath)

# ╔═╡ 8fe1a81f-be2d-4254-8029-aca2dc07a609
seg1 = niread(volseg1);

# ╔═╡ 6aea2d01-b4d5-442a-aafa-fec218cc4123
seg2 = niread(volseg2);

# ╔═╡ 644c1446-4206-4857-8aa3-ca87e0a6ff99
volraw = DICOMUtils.load_dcm_array(vol);

# ╔═╡ d729e367-61e4-4ade-b898-d66631443236
seg1raw = copy(seg1.raw);

# ╔═╡ f883ae8f-5761-471e-b80f-ec55a1fca07a
seg2raw = copy(seg2.raw);

# ╔═╡ 7d4ddb3d-0d5e-433c-a8f1-6ff77cdfd1c5
md"""
### Without Reorientation
"""

# ╔═╡ 5d6edfe9-9d47-4c8b-ac6b-ebc35e3b9915
@bind a Slider(1:size(volraw)[3], default=155, show_value=true)

# ╔═╡ b198dedc-0e5e-4c69-a755-387de3997d09
Plots.heatmap(volraw[:,:,a], c=:grays)

# ╔═╡ 4d49321a-7dbe-428b-91c4-c9561055ce0b
Plots.heatmap(seg1raw[:,:,a], c=:grays)

# ╔═╡ 55744a6c-a224-42cc-b590-bc7ec4614bc7
md"""
### Superimpose
"""

# ╔═╡ d4bc98e5-e92b-4b24-ad36-081b9c0d8507
md"""
## With Reorientation
"""

# ╔═╡ 670ec570-8a19-479d-a619-873c8fe9457a
volraw_ornt, affvol, new_affvol = DICOMUtils.orientation(volraw, (("L", "A", "S")));

# ╔═╡ 3bb89d52-da0e-42d2-9368-41c576b3ef11
seg1raw_ornt, aff1, new_aff1 = DICOMUtils.orientation(seg1raw, (("L", "A", "S")));

# ╔═╡ 16da0125-d6ad-45d8-ba30-431ff2fdc7ed
Plots.heatmap(volraw_ornt[:,:,a], c=:grays)

# ╔═╡ 8c49c4d7-8864-4a5d-9797-3e2bacb2f514
Plots.heatmap(seg1raw_ornt[:,:,a], c=:grays)

# ╔═╡ bb387fef-0034-467f-8219-232a5942c264
md"""
## Calcium Scoring
"""

# ╔═╡ 64599dbe-4f10-476b-921b-63872d7482d1
begin
	ii1 = [1 2 0 1
		   1 5 0 1]
	ii2 = [1 2 0 1
	   	   1 5 0 1]
	img = cat(ii1, ii2, dims=3)
end

# ╔═╡ 228e1866-5d50-41a6-b5f3-cce2d49cdfee
spc = [0.625, 0.625, 0.5]

# ╔═╡ 1252701e-6812-49e4-9282-7c51dd92f8cc
begin
	mm1 = [1 0 0 1
		   1 1 0 1]
	mm2 = [1 0 0 1
	   	   1 1 0 1]
	msk = cat(mm1, mm2, dims=3)
end

# ╔═╡ 288e776b-04ae-479f-8ae2-b3232037fe87
function compute_calcium_scores(image, spacing, mask, min_vol=nothing, max_vol=nothing)
	binary_mask = mask
	voxel_volume = spacing[1] * spacing[2] * spacing[3]
	
	agatston_score = 0
	calcium_volume = 0
	
	# Find individual lesions (in 3D) so that we can discard too small or too large lesions
	lesion_map = Images.label_components(binary_mask)
	n_lesions = length(unique(lesion_map))
	
	for lesion = 1:(n_lesions - 1)
		lesion_mask = map(x -> x == lesion, lesion_map)
		
		# Ignore too small or too large lesions
		lesion_volume = count(lesion_mask) * voxel_volume
		if ((min_vol != nothing) && (lesion_volume < min_vol))
			continue
		end
		if ((max_vol != nothing) && (lesion_volume > max_vol))
			continue
		end
		
		calcium_volume += lesion_volume
		
		# Calculate Agatston score for this lesion
		slices = sort(unique(lesion_mask)) .+ 1
		for z = 1:slices[end]
			fragment_mask = lesion_mask[z,:,:]
			n_pixels = count(fragment_mask)
			maximum_intensity = maximum(image[z,:,:][fragment_mask])
            if maximum_intensity < 200
                coefficient = 1
			elseif maximum_intensity < 300
				coefficient = 2
			elseif maximum_intensity < 400
				coefficient = 3
			else
				coefficient = 4
			end
			agatston_score += coefficient * n_pixels
		end
	end
	agatston_score *= spacing[1] / 3.0 * spacing[2] * spacing[3]
	return agatston_score, calcium_volume
end

# ╔═╡ ea475810-4853-4f6f-97bd-d3faf1b646d6
compute_calcium_scores(img, spc, msk)

# ╔═╡ 165084b0-4ebe-4c28-bdae-24278d7497ae
compute_calcium_scores(img, spc, msk, 0)

# ╔═╡ Cell order:
# ╠═47d5b09c-c213-11eb-2472-7b23859c3caa
# ╠═e991c00b-974b-442c-93c6-55a0b12b27ed
# ╟─718b9504-f15f-41d1-b0ea-3346f17b4cca
# ╠═a45f05b9-08d6-4277-9c9c-205f60b43391
# ╠═5e69faec-3f89-4f65-aab1-9429ad1e141b
# ╠═0399cca4-b251-451e-bf8c-96ccd207a1df
# ╟─20b7a9f9-4308-427f-995f-1f29f7f6df23
# ╠═e36c897d-07ae-4f61-8a5d-6dbdc651431a
# ╠═8fe1a81f-be2d-4254-8029-aca2dc07a609
# ╠═6aea2d01-b4d5-442a-aafa-fec218cc4123
# ╠═644c1446-4206-4857-8aa3-ca87e0a6ff99
# ╠═d729e367-61e4-4ade-b898-d66631443236
# ╠═f883ae8f-5761-471e-b80f-ec55a1fca07a
# ╟─7d4ddb3d-0d5e-433c-a8f1-6ff77cdfd1c5
# ╠═5d6edfe9-9d47-4c8b-ac6b-ebc35e3b9915
# ╠═b198dedc-0e5e-4c69-a755-387de3997d09
# ╠═4d49321a-7dbe-428b-91c4-c9561055ce0b
# ╟─55744a6c-a224-42cc-b590-bc7ec4614bc7
# ╟─d4bc98e5-e92b-4b24-ad36-081b9c0d8507
# ╠═670ec570-8a19-479d-a619-873c8fe9457a
# ╠═3bb89d52-da0e-42d2-9368-41c576b3ef11
# ╠═16da0125-d6ad-45d8-ba30-431ff2fdc7ed
# ╠═8c49c4d7-8864-4a5d-9797-3e2bacb2f514
# ╟─bb387fef-0034-467f-8219-232a5942c264
# ╠═64599dbe-4f10-476b-921b-63872d7482d1
# ╠═228e1866-5d50-41a6-b5f3-cce2d49cdfee
# ╠═1252701e-6812-49e4-9282-7c51dd92f8cc
# ╠═288e776b-04ae-479f-8ae2-b3232037fe87
# ╠═ea475810-4853-4f6f-97bd-d3faf1b646d6
# ╠═165084b0-4ebe-4c28-bdae-24278d7497ae
