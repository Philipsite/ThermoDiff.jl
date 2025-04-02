"""
    apparent_G(P, T, ΔfH°, S°, V°, CPbb84_k1, CPbb84_k4, CPbb84_k3, CPbb84_k8, Vtwq_v1, Vtwq_v2, Vtwq_v3, Vtwq_v4)
Calculate the apparent Gibbs free energy (G) using the Berman and Brown (1983) method.

    Parameters
    ----------
        P : AbstractFloat
    Pressure in bar.            
        T : AbstractFloat
    Temperature in Kelvin.
        ΔfH° : AbstractFloat
    Standard enthalpy of formation in J/mol.
        S° : AbstractFloat
    Standard entropy in J/(mol*K).
        V° : AbstractFloat
    Molar volume in J/(mol*bar).
        CPbb84_k1, CPbb84_k4, CPbb84_k3, CPbb84_k8 : AbstractFloat
    Coefficients for the Cp function from Berman and Brown (1984).
        Vtwq_v1, Vtwq_v2, Vtwq_v3, Vtwq_v4 : AbstractFloat
    Coefficients for the simple V function in TWQ from Theriak Domino Guide.

    Returns
    -------
        ∆aG_PT : AbstractFloat
    Apparent Gibbs free energy in J/mol.

    Notes
    -----   
    Benchmarked against theriak using the JUN92d database by calculating the Gibbs free energy of KYANITE at multiple temperatures and pressures.

"""
function apparent_G(P           :: Union{AbstractFloat, AbstractArray},
                    T           :: Union{AbstractFloat, AbstractArray},
                    ΔfH°        :: AbstractFloat,
                    S°          :: AbstractFloat,
                    V°          :: AbstractFloat,
                    CPbb84_k1   :: AbstractFloat,
                    CPbb84_k4   :: AbstractFloat,
                    CPbb84_k3   :: AbstractFloat,
                    CPbb84_k8   :: AbstractFloat,
                    Vtwq_v1     :: AbstractFloat,
                    Vtwq_v2     :: AbstractFloat,
                    Vtwq_v3     :: AbstractFloat,
                    Vtwq_v4     :: AbstractFloat)

    # Integrate Cp function
    ∫Cp =   ∫Cp_BermanBrown1983(T, CPbb84_k1, CPbb84_k4, CPbb84_k3, CPbb84_k8)
    ∫CpT =  ∫CpT_BermanBrown1983(T, CPbb84_k1, CPbb84_k4, CPbb84_k3, CPbb84_k8)

    # Integrate V function
    ∫V = ∫V_Twq(T, P, V°, Vtwq_v1, Vtwq_v2, Vtwq_v3, Vtwq_v4)

    # Print intermediate results for debugging
    # println("∫Cp: ", ∫Cp)
    # println("∫CpT: ", ∫CpT)
    # println("∫V: ", ∫V)

    # Calculate the apparent Gibbs free energy
    ∆aG_PT = ΔfH° .+ ∫Cp .- T .* S° .- T .* ∫CpT .+ ∫V
    return ∆aG_PT
end


function ∫Cp_BermanBrown1983(T:: Union{AbstractFloat, AbstractArray},
                             k1:: AbstractFloat,
                             k4:: AbstractFloat,
                             k3:: AbstractFloat,
                             k8:: AbstractFloat)
    # Constants
    T° = 298.15             # Reference temperature in Kelvin (25 °C)

    # Integrate the Cp function
    return k1 .* (T .- T°) .- k3 .* (1 ./T .- 1/T°) .+ 2 * k4 .* (sqrt.(T) .- sqrt(T°)) .- k8/2 .* (T.^-2 .- T°^-2)
end


function ∫CpT_BermanBrown1983(T:: Union{AbstractFloat, AbstractArray},
                              k1:: AbstractFloat,
                              k4:: AbstractFloat,
                              k3:: AbstractFloat,
                              k8:: AbstractFloat)
    # Constants
    T° = 298.15             # Reference temperature in Kelvin (25 °C)

    # Integrate the Cp function
    return k1 .* log.(T./T°) .- k3/2 .* (T.^-2 .- T°^-2) .- 2 .* k4 .* (1 ./sqrt.(T) .- 1/sqrt(T°)) .- k8/3 .* (T.^-3 .- T°^-3)
end


function ∫V_Twq(T:: Union{AbstractFloat, AbstractArray},
                P:: Union{AbstractFloat, AbstractArray},
                V°:: AbstractFloat,
                v1:: AbstractFloat,
                v2:: AbstractFloat,
                v3:: AbstractFloat,
                v4:: AbstractFloat)
    
    # Constants
    T° = 298.15            # Reference temperature in Kelvin (25 °C)
    P° = 1.0               # Reference pressure in bar

    # convert v1, v2, v3, v4 (see Theriak Domino Guide, p.53)
    v1 = v1 * 1e-5
    v2 = v2 * 1e-5
    v3 = v3 * 1e-5
    v4 = v4 * 1e-8

    # precalculation
    Vta = v1 * V°
    Vtb = v2 * V°
    Vpa = v3 * V°
    Vpb = v4 * V°

    # Integrate the V function
    return (V° .+ Vta .* (T .- T°) .+ Vtb .* (T .- T°).^2) .* (P .- P°) .+
            Vpa .* (P.^2 ./ 2 .- P .* P° .+ P°.^2 ./ 2) .+ Vpb .* (P.^3 ./ 3 .- P.^2 .* P° .+ P .* P°.^2 .* P°.^3 ./ 3)
end
