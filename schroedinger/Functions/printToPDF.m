function printToPDF(myHandle,name)
% Input:
%   myHandle    - The handle to the figure which is to be saved
%   name        - A string containing the file name

%     screen_size = get(0, 'ScreenSize');
%     origSize = get(myHandle, 'Position'); % grab original on screen size
%     set(myHandle, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to scren size
%     set(myHandle,'PaperPositionMode','auto') %set paper pos for printing
%     print(myHandle,name,'-dpdf','-r0') % save figure
%     set(myHandle,'Position', origSize) %set back to original dimensions

    set(myHandle,'Units','Inches');
    pos = get(myHandle,'Position');
    set(myHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(myHandle,name,'-dpdf','-r800')
    pause(0.1)
    % Move to plots folder. If folder does not exist, create it.
    if ~exist('Plots', 'dir')
        mkdir('Plots');
    end
    movefile([name '.pdf'],'Plots')
end