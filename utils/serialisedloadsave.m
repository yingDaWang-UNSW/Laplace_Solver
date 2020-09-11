function [data]=serialisedloadsave(fileName,data,saveFlag,serialiseFlag)
    if saveFlag
        if ~serialiseFlag
        save(fileName,'data','-v7.3');
        else
            disp('Serialising data')
            data=hlp_serialize(data);
            savefast(fileName,'data');
            data=[];
        end
    else
        if ~serialiseFlag
        data=load(fileName);
        else
            load(fileName,'data');
            disp('Deserialising data')
            data=hlp_deserialize(data);
        end
    end
end