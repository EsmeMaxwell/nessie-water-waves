function fn_ww__setup__prepare_mp( st_Dn, st_p )
%FN_WW__SETUP__PREPARE_MP: Prepare Advanpix MP (may need to do per-process)
%
%   FN_WW__SETUP__PREPARE_MP()
%

if ( st_p.bp_mp )
    assert( ~isfield( st_p, 'ip_mp_digits' ) || ( isfield( st_Dn, 'mp_digits' ) && st_Dn.mp_digits == st_p.ip_mp_digits ), 'Inconsistent setting of mp.Digits' )
    
    if ( mp.Digits() ~= st_Dn.mp_digits )
        mp.Digits( st_Dn.mp_digits );
    end    
    
    if ( 0 ~= st_p.ip_mp_threads )
        mp.NumberOfThreads( st_p.ip_mp_threads );
    end
end


end