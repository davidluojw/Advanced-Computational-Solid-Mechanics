function [val] = f(xx, yy)

dist = sqrt( (xx-0.5).^2 + (yy-0.5).^2 );
if dist < 0.05
  val = 1.0;
else
  val = 0.0;
end