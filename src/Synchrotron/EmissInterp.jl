"""
    Equations from Padovani et. 2021 (https://doi.org/10.1051/0004-6361/202140799) for the interpolation code 
"""

nu_c(E,BField) = PRE_NU_C * BField * (E / ELECTRON_ENERGY_AT_REST_eV)^2

freq_norm(E,nu,BField) = nu / (nu_c(E,BField))

F_integrand(x) = besselk(5/3,x)

F(x) = x * quadgk(F_integrand,x,Inf, rtol=1.e-3)[1]

G(x) = x * (besselk(2/3,x))

relativistic_ve(E) = C*(sqrt(1-(1/(1+E/ELECTRON_ENERGY_AT_REST_eV)^2)))

j_e(E,a=-1.3,b=1.9) = J_0 * E^a / (E+E_0)^b #Padovani 2021

je_ve_ratio(E) = j_e(E) / relativistic_ve(E)

power_par(E,nu,BField) = PRE_P * BField * (F(freq_norm(E,nu,BField)) - G(freq_norm(E,nu,BField)))
power_perp(E,nu,BField)= PRE_P * BField * (F(freq_norm(E,nu,BField)) + G(freq_norm(E,nu,BField)))

par_integrand(E,nu,BField) = je_ve_ratio(E) * power_par(E,nu,BField)
perp_integrand(E,nu,BField) = je_ve_ratio(E) * power_perp(E,nu,BField)

function par_emissivity(nu_MHz,BField_microG)
    nu = 1e6 * nu_MHz
    BField = 1e-6 * BField_microG
    return quadgk(x -> par_integrand(x,nu,BField),ELECTRON_ENERGY_AT_REST_eV,1e10)[1]
end

function perp_emissivity(nu_MHz,BField_microG)
    nu = 1e6 * nu_MHz
    BField = 1e-6 * BField_microG
    return quadgk(x -> perp_integrand(x,nu,BField),ELECTRON_ENERGY_AT_REST_eV,1e10)[1]
end

"""
    EmissInterp(BArray::AbstractArray, nuArray::AbstractArray)

Calculate the emissivity for a range of magnetic fields and frequencies, and write the results to a file.

# Arguments
- `BArray::AbstractArray`: Array of magnetic field strengths in microGauss.
- `nuArray::AbstractArray`: Array of frequencies in MHz.

# Returns
- `Nothing`: This function does not return a value but writes the results to a file named "emissivite_interp_LOFAR.dat".

# Description
This function calculates the parallel and perpendicular emissivities for each combination of magnetic field strengths and frequencies provided in `BArray` and `nuArray`. The results are written to a file named "emissivity.dat" in a tab-separated format with columns for magnetic field strength (`B`), frequency (`nu`), parallel emissivity (`e_para`), and perpendicular emissivity (`e_perp`).

# Example
```julia
BArray = [1.0, 2.0, 3.0]
nuArray = [100, 200, 300]
EmissInterp(BArray, nuArray)
"""

function EmissInterp(BArray::AbstractArray,nuArray::AbstractArray)
    open("emissivity.dat", "w") do f
      write(f, "B\tnu\te_para\te_perp\n")    
      for nui in nuArray   
        for Bi in BArray
            @inbounds eps_para = par_emissivity(nui,Bi)
            @inbounds eps_perp = perp_emissivity(nui,Bi)
            to_print = "$Bi\t$nui\t$eps_para\t$eps_perp\n"
            write(f, to_print)
        end
      end
    end
end