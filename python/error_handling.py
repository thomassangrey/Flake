def WARN_lswt_too_small(QL_mask, unc_mask, lswt_mask, QL, unc_lswt, lswt, Q_mask, t_idx):
    print(f'WARNING": for t_idx = {t_idx}, lswt is too small.')
    print(f'lswt.size must be > Q_mask.size')
    print(f'    QL_mask.shape = {QL_mask.shape}; unc_mask.shape = {unc_mask.shape};\n')        
    print(f'    lswt_mask.shape = {lswt_mask.shape}\n')
    print(f'    QL.shape = {QL.shape}; unc_lswt.shape = {unc_lswt.shape}; \n')
    print(f'    lswt.shape = {lswt.shape}; Q_mask.shape = {Q_mask.shape})\n')
    err = True
    return err

def WARN_Q_mask_no_pixels(not_warned):
    if not_warned:
        not_warned = False
        print('WARNING: Quality mask is entirely True (not enough quality pixels). Returning Nan.\n')
    return not_warned