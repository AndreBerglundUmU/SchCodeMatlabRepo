function printToPDF(myHandle,name)
% printToPDF    - Save the figure to PDF format in the folder "Plots".
% Syntax: printToPDF(myHandle,name)
%
% Input:
% myHandle  - The handle to the figure which is to be saved
% name      - A string containing the file name
%
% Output:
% A file 'name.pdf' in the folder 'Plots' with quality level r800.
%
% Non-standard dependencies: None.
% See also: Any accompanying script for example usage.
%           printToEPS.m

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