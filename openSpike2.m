function info = openSpike2(path)

s_FileID    = fopen(path);


if s_FileID < 0
    error('[openSpike2] - No FileID')
end

s_Channels  = 32;
s_SampScale	= 1e-6;
info = repmat(struct('data',[],'header',[]),s_Channels,1);


for s_cCh = 1:s_Channels
    
    [v_Data,st_Header]	= SONGetChannel(s_FileID,s_cCh);
           
    if ~isempty(v_Data)
            if s_cCh ==1 || s_cCh == 2
            [v_Data,st_Header]	= SONADCToDouble(v_Data,st_Header);
            v_Data              = v_Data(:)';
             st_Header.Sampling  = 1/(st_Header.sampleinterval*s_SampScale);
            end
         info(s_cCh).data    = v_Data;
         info(s_cCh).header  = st_Header;
    else continue
    end
end
fclose(s_FileID);
end