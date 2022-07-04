%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% permited_it
% Returns a boolean if the proposed iteration is allowed or not.
% Customizable for each dynamic model.
%
% INPUTS:
% it            Structure with all iterations. Details of this structure
%               can be found in HIPPO.m.
% last_results  Structure with the last iteration. Details of this
%               structure can be found in iteration.m.
% 
% OUTPUT:
% answer        If the iteration is permited or not.
%
% NOTE: The new iteration is in it.new, which is a vector with ones if the
% parameter is fixed, and zeros if the parameter is still free.
%
% Benjamín J. Sánchez
% Last Update: 2014-07-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function answer = permited_it(it,last_results)

%Ignore solutions with unexpected sensitivity outputs:
if max(max(last_results.Ms)) >= 7
    answer = false;
       
%Problem-specific constraints (solid-substrate fermentation):

% %[th2 th3 th4 th5 th12 th16 th17]
% [pilipmax picitmax Yxs ks2 k2 rhonmax Q0]

% ds = Sin*(fsin/V)-((mu/Yxs)+(pilip/Ylips)+(picit/Ycits))*Xf-(S/V)*Fin;
elseif it.new(1) + it.new(2) + it.new(3) + it.new(4) + it.new(5) + it.new(7) == 6
    answer = false;                                        
                                                           
%dn = Nin*(fnin/V)-rhon*Xf-(N/V)*Fin;    
elseif it.new(6) == 1
    answer = false;                                        

%dxf = mu*Xf-(Xf/V)*Fin;
elseif it.new(7) == 1
    answer = false;
% %dxl = (pilip+gamma*mu)*Xf-(Xl/V)*Fin;
% elseif it.new(1) + it.new(17) + it.new(7) == 3
%     answer = false;

%dcit = picit*Xf-(CIT/V)*Fin;
elseif it.new(2) + it.new(5) == 1
    answer = false;
%dqn = rhon-mu*qn;
elseif it.new(6) + it.new(7) == 2
    answer = false;

% [pilipmax picitmax Yxs ks2 k2 rhonmax Q0]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DO NOT MODIFY THIS SECTION
else
    U      = evalin('base','U');
    answer = true;
    %Ignore solution if last_results has more problems than his father:
    if sum(it.last) ~= 0
        problems_child = sum(last_results.ktofix);
        for i = 1:length(it.codes)
            R    =  it.last - it.codes{i,1};
            diff = sum(R.^2);
            if diff ==0
                pos_child = i;
            end
        end
        for i = 1:length(it.tree(:,1))
            if it.tree(i,pos_child) == 1
                problems_father = sum(it.codes{i,2}.ktofix);
                if problems_child > U*problems_father
                    answer = false;
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%