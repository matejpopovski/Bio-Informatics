def lacI_bound(lactose, lacI):
    return not lactose and lacI

def CAP_bound(glucose, CAP):
    return not glucose and CAP

def lacZ(is_lacI_bound, is_CAP_bound):
    if is_lacI_bound == False and is_CAP_bound == False:
        return "low"
    elif is_lacI_bound == False and is_CAP_bound == True:
        return "high"
    else:
        return "absent"
def lacZ_full_circuit(lactose, lacI, glucose, CAP):
    return lacZ(lacI_bound(lactose, lacI), CAP_bound(glucose, CAP))