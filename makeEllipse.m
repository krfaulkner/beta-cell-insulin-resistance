% function: makeEllipse makes a filled ellipse

function el=makeEllipse(a,b,x,y,borderColor,fillColor)
    t=linspace(0,2*pi,100);
    xs=b*cos(t)+x;
    ys=a*sin(t)+y;
    el=patch(xs,ys, fillColor,'EdgeColor', borderColor,'FaceAlpha',0.2);
end