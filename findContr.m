function [isContr, contr] = findContr(x1,y1,x2,y2,Dx,Dy,X,Y)

X_pt = [X-Dx/2; X+Dx/2];
Y_pt = [Y-Dy/2; Y+Dy/2];
list = [X_pt(1) Y_pt(1); X_pt(2) Y_pt(1); X_pt(2) Y_pt(2); ...
        X_pt(1) Y_pt(2); X_pt(1) Y_pt(1)];

% finding if it intersects
% calculating where max and min x and y fall
tx = (X_pt-x1)./(x2-x1);
ty = (Y_pt-y1)./(y2-y1);
% if there is an overlap in segment tx and ty then there is an intersection
txy=[tx;ty];
if (~(all(txy<0) || all(txy>1)))
    tx1=min(tx); tx2=max(tx);
    ty1=min(ty); ty2=max(ty);
    if((ty1>tx2) || (tx1>ty2)>0)
        isContr=0; contr=0;
        return;
    else
        isContr=1;
    end
else
    isContr=0; contr=0;
    return;
end

% case one if any of the end points are inside
% case two if any of the sides intersects with the segment
% hit_bnd=0;
if isContr==1
    if (x1>X_pt(1) && x1<X_pt(2) && y1>Y_pt(1) && y1<Y_pt(2)) % starting point is inside
        if (x2>X_pt(1) && x2<X_pt(2) && y2>Y_pt(1) && y2<Y_pt(2))
            % both points are inside
            %Contr_typ=1;
%             isContr=1;
            contr = sqrt((x1-x2)^2 + (y1-y2)^2);
        else
            % end point is outside but start point is in
            Contr_typ=0; % resetting the value
            for ii=1:4
                x3 = list(ii,1); x4=list(ii+1,1);
                y3 = list(ii,2); y4=list(ii+1,2);
                ax=x2-x1; ay=y2-y1;
                bx=x4-x3; by=y4-y3;
                det=ax*by-ay*bx;
                
                if(det~=0)
                    t = ((y1-y3)*(x4-x3) - (x1-x3)*(y4-y3))/det;
                    s=-((y3-y1)*(x2-x1) - (x3-x1)*(y2-y1))/det;
                    
                    if(t>=0 && t<=1 && s>=0 && s<=1)
                        if(det>0)
                            Contr_typ=1; %hit_bnd=ii;
                            %                     else
                            %                         isContr=-1; hit_bnd=ii;
                        end
                    end
                end
                if Contr_typ==1
                    break;
                end
            end
            x_int = x1 + (x2-x1)*t;
            y_int = y1 + (y2-y1)*t;
            % distance between starting point and intersection point
            isContr=1;
            contr = sqrt((x_int-x1)^2 + (y_int-y1)^2);
        end
    else
        if(x2>X_pt(1) && x2<X_pt(2) && y2>Y_pt(1) && y2<Y_pt(2))
            % start point is out but end point is in
            % now we are searching for negative determinant
            Contr_typ=0; % resetting the value
            for ii=1:4
                x3 = list(ii,1); x4=list(ii+1,1);
                y3 = list(ii,2); y4=list(ii+1,2);
                ax=x2-x1; ay=y2-y1;
                bx=x4-x3; by=y4-y3;
                det=ax*by-ay*bx;
                
                if(det~=0)
                    t = ((y1-y3)*(x4-x3) - (x1-x3)*(y4-y3))/det;
                    s=-((y3-y1)*(x2-x1) - (x3-x1)*(y2-y1))/det;
                    
                    if(t>=0 && t<=1 && s>=0 && s<=1)
                        if(det<0)
                            Contr_typ=-1; %hit_bnd=ii;
                        end
                    end
                end
                if Contr_typ==-1
                    break;
                end
            end
            x_int = x1 + (x2-x1)*t;
            y_int = y1 + (y2-y1)*t;
            % distance between intersection point and final point
%             isContr=1;
            contr = sqrt((x2-x_int)^2 + (y2-y_int)^2);
            
        else
            % both points are out
            % here we need two points one with positive determinant another
            % with negative
            % first negative one to get starting point
            Contr_typ=0; % resetting the value
            for ii=1:4
                x3 = list(ii,1); x4=list(ii+1,1);
                y3 = list(ii,2); y4=list(ii+1,2);
                ax=x2-x1; ay=y2-y1;
                bx=x4-x3; by=y4-y3;
                det=ax*by-ay*bx;
                
                if(det~=0)
                    t = ((y1-y3)*(x4-x3) - (x1-x3)*(y4-y3))/det;
                    s=-((y3-y1)*(x2-x1) - (x3-x1)*(y2-y1))/det;
                    
                    if(t>=0 && t<=1 && s>=0 && s<=1)
                        if(det<0)
                            Contr_typ=-1; %hit_bnd=ii;
                        end
                    end
                end
                if Contr_typ==-1
                    break;
                end
            end
            x_int1 = x1 + (x2-x1)*t;
            y_int1 = y1 + (y2-y1)*t;
            
            % now the positive one to get the final point
            Contr_typ=0; % resetting the value
            for ii=1:4
                x3 = list(ii,1); x4=list(ii+1,1);
                y3 = list(ii,2); y4=list(ii+1,2);
                ax=x2-x1; ay=y2-y1;
                bx=x4-x3; by=y4-y3;
                det=ax*by-ay*bx;
                
                if(det~=0)
                    t = ((y1-y3)*(x4-x3) - (x1-x3)*(y4-y3))/det;
                    s=-((y3-y1)*(x2-x1) - (x3-x1)*(y2-y1))/det;
                    
                    if(t>=0 && t<=1 && s>=0 && s<=1)
                        if(det>0)
                            Contr_typ=1; %hit_bnd=ii;
                        end
                    end
                end
                if Contr_typ==1
                    break;
                end
            end
            x_int2 = x1 + (x2-x1)*t;
            y_int2 = y1 + (y2-y1)*t;
            
%             isContr = 1;
            contr = sqrt((x_int1-x_int2)^2 + (y_int1-y_int2)^2);
        end
    end
  
end
