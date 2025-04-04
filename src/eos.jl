"""
    apparent_G(P, T, ΔfH°, S°, V°, CPbb84_k1, CPbb84_k4, CPbb84_k3, CPbb84_k8, Vtwq_v1, Vtwq_v2, Vtwq_v3, Vtwq_v4)
Calculate the apparent Gibbs free energy (G) using the Berman and Brown (1983) method.

# Parameters
    P : Pressure in bar.            
    T : Temperature in Kelvin.
    ΔfH° : Standard enthalpy of formation in J/mol.
    S° : Standard entropy in J/(mol*K).
    V° : Molar volume in J/(mol*bar).
    CPbb84_k1, CPbb84_k4, CPbb84_k3, CPbb84_k8 : Coefficients for the Cp function from Berman and Brown (1984).
    Vtwq_v1, Vtwq_v2, Vtwq_v3, Vtwq_v4 : Coefficients for the simple V function in TWQ from Theriak Domino Guide.

# Returns
    ∆aG_PT : Apparent Gibbs free energy in J/mol.

# Notes 
Benchmarked against theriak using the JUN92d database by calculating the Gibbs free energy of KYANITE at multiple temperatures and pressures.

"""
function apparent_G(P           :: Union{Number, AbstractArray},
                    T           :: Union{Number, AbstractArray},
                    ΔfH°        :: Number,
                    S°          :: Number,
                    V°          :: Number,
                    CPbb84_k1   :: Number,
                    CPbb84_k4   :: Number,
                    CPbb84_k3   :: Number,
                    CPbb84_k8   :: Number,
                    Vtwq_v1     :: Number,
                    Vtwq_v2     :: Number,
                    Vtwq_v3     :: Number,
                    Vtwq_v4     :: Number)

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


function ∫Cp_BermanBrown1983(T:: Union{Number, AbstractArray},
                             k1:: Number,
                             k4:: Number,
                             k3:: Number,
                             k8:: Number)
    # Constants
    T° = 298.15             # Reference temperature in Kelvin (25 °C)

    # Integrate the Cp function
    return k1 .* (T .- T°) .- k3 .* (1 ./T .- 1/T°) .+ 2 * k4 .* (sqrt.(T) .- sqrt(T°)) .- k8/2 .* (T.^-2 .- T°^-2)
end


function ∫CpT_BermanBrown1983(T:: Union{Number, AbstractArray},
                              k1:: Number,
                              k4:: Number,
                              k3:: Number,
                              k8:: Number)
    # Constants
    T° = 298.15             # Reference temperature in Kelvin (25 °C)

    # Integrate the Cp function
    return k1 .* log.(T./T°) .- k3/2 .* (T.^-2 .- T°^-2) .- 2 .* k4 .* (1 ./sqrt.(T) .- 1/sqrt(T°)) .- k8/3 .* (T.^-3 .- T°^-3)
end


function ∫V_Twq(T:: Union{Number, AbstractArray},
                P:: Union{Number, AbstractArray},
                V°:: Union{Number, Any},
                v1:: Number,
                v2:: Number,
                v3:: Number,
                v4:: Number)
    
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
