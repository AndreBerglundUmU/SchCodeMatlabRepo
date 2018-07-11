function printToEPS(myHandle,name)
% printToPDF    - A function which saves the figure to EPS format in the
%                 folder "Plots", using quality level r800.
% Syntax: printToPDF(myHandle,name)
%
% Input:
% myHandle  - The handle to the figure which is to be saved
% name      - A string containing the file name
%
% Output:
% A file 'name.eps' in the folder 'Plots'.
%
% Non-standard dependencies: None.
% See also: printToPDF.m

    set(myHandle,'Units','Inches');
    pos = get(myHandle,'Position');
    set(myHandle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print(myHandle,name,'-depsc','-r800')
    pause(0.1)
    % Move to plots folder. If folder does not exist, create it.
    if ~exist('Plots', 'dir')
        mkdir('Plots');
    end
    movefile([name '.eps'],'Plots')
end