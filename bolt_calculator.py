from math import *
import numpy as np
import pdb
# NOTE: Imperial Units in, lb
# NOTE: it would be helpful to have a picture associated with the calculations, that highlights which part of the frustum that is calculated in the thickness
# NOTE: Refactor and remove global variables if necessary
procedure = 1

Em1 = 2.85e+07
Ewasher = 30e6
Em2 = 2.9e+07
H = 0 #nut thickness

d = 3/8
dw = 1.5*d 
db = 1.312 #diameter of bolt
tw = 0.083
tw2 = 0
tm1 = 0.1345
tm2 = 0.25

At = 0.0775
nominal_diameter = 0.375
minor_diameter = 0.0678
E_bolt = 30e6
ld = 0 #length of unthreaded portion; zero by default
L = 0 #shank length of the bolt
def calculate_washer_plus_material1_thickness(tw, tm1):
    if tw:
        return tw + tm1 #thickness of washer plus thickness of material 1
    return 0
def grip_length(procedure, t2=0, h=0, d=0):
    match procedure:
        case 1:
            return tw + tw2 + tm1 + tm2
        case 2:
            if t2 < d:
                return h + t2/2
            return h + d/2
        case _:
            raise Exception("Please choose a procedure")
def fastener_length(procedure, l=0, H=0, h=0, d=0):
    match procedure:
        case 1:
            print(f"Minimum Shank Length: {l + H} in")
            return l + H
        case 2:
            print(f"Minimum Shank Length: {h + 1.5*d} in")
            return h + 1.5*d
        case _:
            raise Exception("Please choose a procedure")
def upper_frusta(Ewasher, Em1):
    #if E == Ew then thickness tw + tm1
    #else then return a different calculation
    if Ewasher == Em1:
        return  {"name": "washer 1 with material 1","t": tw + tm1, "D": dw, "E": Em1}
    return [
        {"name": "top washer", "t": tw, "D": dw, "E": Ewasher}, 
        {"name": "material 1", "t": tm1, "D": dw + 2*tw*tan(radians(30)), "E": Em1}
    ]
def middle_frusta(g_length, E):
    if tw:
        return {
            "name": "midpoint frustum material",
            "t": g_length/2 - tm1 - tw,
            "D": dw + 2*(tw + tm1) * tan(radians(30)),
            "E": E
        }
    if tm1 != tm2:
        return {
            "name": "midpoint frustum material",
            "t": g_length/2 - tm1,
            "D": dw + 2*tm1 * tan(radians(30)), #even this there is no washer, dw is a close estimate to the diameter of the bolt head
            "E": E
        }
    return { #NOTE: There is no midpoint frustum in this case
        "name": "midpoint frustum material",
        "t": 0,
        "D": 0, 
        "E": 0
    }
    
def lower_frusta(g_length, E, Ewasher=0):
    if tw2 and g_length >= (tw + tw2 + tm1 + tm2): #if the g length extends through the washer
        if E == Ewasher or Ewasher == 0:
            return {"name": "material 2", "t": g_length/2, "D": dw, "E": E}
        return [
            {"name": "material 2","t": grip_length/2 - tw2, "D": dw + 2*tw2*tan(radians(30)), "E": E}, 
            {"name": "washer 2", "t": tw2, "D": dw, "E": Ewasher}
        ]
        
    return {"name": "material 2","t": g_length/2,"D": dw, "E": E} #NOTE: when dw == 0, then dw = db
def flatten_array(arr):
    out = np.array([])
    for i in frustum:
        out = np.append(out, i).flatten()
    return out
def stiffness(E, t, D, d):
    if t == 0:
        print(f"NOTE: t is zero {t}")
        return 0 
    #case x:
    numerator = 0.5774*pi*d*E

    ln_numerator = (1.155*t + D - d) * (D + d)
    ln_denominator = (1.155*t + D + d) * (D - d)

    ln_qoutient = log(ln_numerator/ln_denominator)
    
    return numerator/ln_qoutient
def threaded_length(L, d): #-> LT
    if L <= 6: #6inches
        return 2*d + 0.25
    return 2*d+0.5
def length_of_unthreaded_portion_in_grip(L, LT): # -> ld
    return L - LT
def length_of_threaded_in_grip(l, ld): # -> lt
    return l - ld
def area_of_unthreaded_portion(d):
    return (pi/4) * d**2
def calculate_stiffnesses(frustum):
    ks = []
    for i in frustum:
        e = i["E"]
        t = i["t"]
        D = i["D"]
        
        ks.append(stiffness(e, t, D, d))
    return ks
def convert_equivalent_stiffness(ks):
    for i, v in enumerate(ks):
        if v == 0:
            continue
        ks[i] = 1/v
    return 1/sum(ks)
def calculate_bolt_stiffness(g_length, procedure, E_bolt, H, h, d, At, L=0, ld=0):
    l = g_length

    if not L:
        L = fastener_length(procedure, g_length, H, h, d) #inch series, #NOTE: this is supposed to be rounded up to a standard fractional length of a bolt. so if 0.78, round up to 7/8 instead of 3/4
    LT = threaded_length(L, d)
    if LT > L:
        L = LT
    if not ld:
        ld = length_of_unthreaded_portion_in_grip(L, LT)
    
    lt = g_length - ld
    Ad = area_of_unthreaded_portion(d)
    
    if ld > 0:
        return (Ad*At*E_bolt) / (Ad * lt + At*ld)
    
    return At * E_bolt / lt
def estimated_stiffness():
    #using steel
    A = 0.7815
    B = 0.62873
    km = (30e6 * d * A * exp(B * d / g_length))/(10**6)
    return km

if __name__ == "__main__":
    h = calculate_washer_plus_material1_thickness(tw,tm1)
    g_length = grip_length(procedure, tm2, h, d)
    if dw == 0:
        dw = db 
   
    kb = calculate_bolt_stiffness(g_length, procedure, E_bolt, H, h, d, At, L)

    frustum = [upper_frusta(Ewasher, Em1), middle_frusta(g_length, Em2), lower_frusta(g_length, Em2)]
    frustum = flatten_array(frustum)
    
    ks = calculate_stiffnesses(frustum)
    km = convert_equivalent_stiffness(ks)
    
    SF = 1
    #P = 73 * 9 * SF #the axial force from moment forces is 73lb, the force multiplier it 9
    P = 73*9*1.5
    C = kb/(kb+km)
    Pb = C * P
    
    Sp = 120e3
    Sut = 150e3
    Se = 23.2e3

    N = 3
    Fi = 0.75*At*Sp
    Fb = Pb + Fi
    sigma_b = Fb / At

    load_factor = (Sp*At - Fi) / (C *(P/N))
    yield_factor = (Sp*At) / (C*(P/N) + Fi)
    joint_seperation_factor = Fi / ((P/N)*(1-C))    
    
    K = 0.2

    Pmin = 0
    Pmax = P

    Fbmin = C*Pmin + Fi
    Fbmax = C*Pmax + Fi
    
    sigma_i = Fi/At
    sigma_a = C*(Pmax + Pmin)/(2 * At)
    sigma_m = sigma_a + Fi/At
    
    fatigue_goodman_SF = Se*(Sut - sigma_i)/ (sigma_a * (Sut + Se))
    fatigue_yield_SF = Sp/(sigma_a + sigma_m)
    

    bolt_data = {
        "Joint Stiffness": f"{round(km*1e-6,2)} Mlbf/in",
        "Bolt Stiffness": f"{round(kb*1e-6,2)} Mlbf/in",
        "Bolt Pretension at 75% Proof Strength": f"{round(Fi)} lb",
        "Minimum Torque": f"{round(Pb * K * d/12, 2)} lb*ft",
        "Recommended Torque": f"{round(K*Fi*d/12)} lb*ft",
        "External Force": f"{round(P)} lb",
        "Stiffness Coefficient (C)": f"{round(C, 2)}",
        "Bolt External Force With Stiffness Modification": f"{round(Pb)} lb",
        "Bolt Total Tension": f"{round(Fb)} lb",
        
        "Bolt Stress": f"{round(sigma_b*1e-3)} ksi",
        "Alternating Stress": f"{round(sigma_a*1e-3)} ksi",
        "Midrange Stresses": f"{round(sigma_m*1e-3)} ksi",

        "Safety Factor": f"{SF}",
        "Load Factor": f"{round(load_factor,2)}",
        "Yield Factor": f"{round(yield_factor,2)}",
        "Joint Separation Factor": f"{round(joint_seperation_factor, 2)}",
        "Fatigue factor of Safety (Goodman Criteria) (Must be greater than 1 to avoid failure)": f"{round(fatigue_goodman_SF, 2)}",
        "Fatigue Yield Factor of Safety": f"{round(fatigue_yield_SF, 2)}" #This determines if it yields at all  under the fatigue loading
    }

    for k,v in bolt_data.items():
        print(f"{k}: {v}")
    print()
    
    