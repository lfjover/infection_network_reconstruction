function setfigure(xSize,ySize,xpos,ypos)
% creates a figure with size X,Y in position xpos ypos and prints the same
% size (units are cm).

set(gcf, 'Units','centimeters', 'Position',[xpos ypos xSize ySize])
set(gcf, 'PaperUnits','centimeters')
set(gcf, 'PaperSize',[xSize ySize])
set(gcf, 'PaperPosition',[0 0 xSize ySize])