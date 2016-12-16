% checks whether a given upstate has action potentials and if asked removes
% them
% inputs :
% - v_Vm : intracellular signal 
% - s_sampleinterval : sampling interval of the signal
% - b_removeAP : if true, action potentials will be removed if found
% outputs :
% - bool : true if length of the upstate signal is more than 300 ms and
% there are no action potentials
% - v_VmPost : signal without action potentials if b_removeAp is true,
% input signal otherwise

function [bool,v_VmPost] = validity_test(v_Vm,...
        s_sampleinterval,b_removeAP)
    
    if nargin == 1
            error('[validity_test] - ERROR: Sample interval missing')   
    end
    
    if nargin == 2
        b_removeAP = false;
    end
    
    b_len = length(v_Vm)*s_sampleinterval/1000 > 300;
    v_diff = (v_Vm(2:end) - v_Vm(1:end-1))/360*1000;
    
    if b_removeAP
        
        b_pa = sum(v_diff >10) <1.5;
        
        if b_pa && sum(v_diff >10) >0.5
             v_VmPost = traitementpa(v_Vm,s_sampleinterval); 
        else v_VmPost = v_Vm;
        end
        
    else b_pa = sum(v_diff >10) <0.5;
        v_VmPost = v_Vm;
    end
    
    bool = b_len && b_pa;
    
end