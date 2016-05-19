function retval = isCommandWindowOpen()
%isCommandWindowOpen check whether running in GUI
    jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
    retval = ~isempty(jDesktop.getClient('Command Window'));
end


