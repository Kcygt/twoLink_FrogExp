function [sX, sY,cX,cY] = defineSquare(x_c, y_c, a)

    % Calculate half side length
    half_side = a / 2;

    % Define vertices
    cX = [ x_c - half_side,  x_c - half_side,  x_c + half_side, x_c + half_side, x_c - half_side];
    cY = [ y_c - half_side,  y_c + half_side,  y_c + half_side, y_c - half_side, y_c - half_side];

    mX = [(cX(2) + cX(1))/2, (cX(3) + cX(2))/2, (cX(4) + cX(3))/2, (cX(5) + cX(4))/2];
    mY = [(cY(2) + cY(1))/2, (cY(3) + cY(2))/2, (cY(4) + cY(3))/2, (cY(5) + cY(4))/2];
    
    sX = [cX(1), mX(1),cX(2), mX(2),cX(3), mX(3),cX(4), mX(4),cX(5)];
    sY = [cY(1), mY(1),cY(2), mY(2),cY(3), mY(3),cY(4), mY(4),cY(5)];
    
end
