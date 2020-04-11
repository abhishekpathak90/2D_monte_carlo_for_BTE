function [xout, yout, frac_in, hit_bnd, scat_type] = Check_intersect(x1,y1,x2,y2,Segments)

[no_seg, ~] = size(Segments);
% ref--- http://www.cs.swan.ac.uk/~cssimon/line_intersection.html


% loop to check for intersection with segments
% based on Cormen chapter 33 (page 937)
% intersect = false;
% hit_bnd = 0;
% ii=0;
% while(~intersect)
%     ii=ii+1;
%     x3=Segments(ii,1); y3=Segments(ii,2);
%     x4=Segments(ii,3); y4=Segments(ii,4);
%     
%     intersect = segment_intersect(x1,y1,x2,y2,x3,y3,x4,y4);
% end
% 
% if(intersect)
%     hit_bnd = ii;
% end

for i=1:no_seg
    x3=Segments(i,1); y3=Segments(i,2);
    x4=Segments(i,3); y4=Segments(i,4);
    
    
    A = [x4-x3 x1-x2; y4-y3 y1-y2];
    B = [x1-x3; y1-y3];
    
    if(det(A) ~=0)
        Check_seg = A\B;
        t2 = Check_seg(1);
        t1 = Check_seg(2);
    end
    
%     Denom1 = (x4-x3)*(y1-y2) - (x1-x2)*(y4-y3);
%     Neum1 = (y3-y4)*(x1-x3) + (x4-x3)*(y1-y3);
%     
%     Denom2 = Denom1;
%     Neum2 = (y1-y2)*(x1-x3) + (x2-x1)*(y1-y3);
%     
%     if(Denom1==0)
%         continue;
%     end
%     
%     t1 = Neum1/Denom1;
%     t2 = Neum2/Denom2;
        
    % This excludes the cases where particle falls on one of the corners.
    % That case is still problematic and needs to be handled later
    if(t1>1e-8 && t1<1 && t2>1e-8 && t2<1)
        hit_bnd = i;
        scat_type = Segments(i,5) + 1;
        
        xout = x1 + t1*(x2-x1);
        yout = y1 + t1*(y2-y1);
        
        frac_in = t1;
        break;
    end
end