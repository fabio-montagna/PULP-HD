%This script creates a file "data.h" that can be incloded in STM/PULP
%projects. 

%select data form subject 1 (COMPLETE_1', hdc_model_1.value), 
%subject 2 (COMPLETE_2', hdc_model_2.value), 
%subject 3 (COMPLETE_3', hdc_model_3.value), etc. 
test_set = COMPLETE_1';
AM = hdc_model_1.values;

%you need to downsample test_set to let it fit in memory. 
downsample = 10; 
test_set=test_set(:,1:downsample:end);
IM= iMch.values;
IM = cell2mat(IM');
CIM = chAM.values;
CIM = cell2mat(CIM'); 
AM = cell2mat(AM');
IM = compress_hypervectors(IM);
CIM = compress_hypervectors(CIM);
AM = compress_hypervectors(AM);

ID = fopen ( 'data.txt' , 'w');
fprintf(ID, '#ifndef DATA_H_\n#define DATA_H_\n\n#include <stdio.h>\n#include "init.h"\n\n//EMG envelope signals\n', r, c);

[r,c]=size(test_set);
fprintf(ID, 'const float TEST_SET[%d][%d] = {\n', r, c);
 for i = 1 : r
    fprintf(ID, '{');
    for j = 1 : c
        if j == c
            fprintf(ID, '%f',(test_set(i, j)));
        else 
            fprintf(ID, '%f, ',(test_set(i, j)));  
        end
    end
    
    fprintf(ID, '},\n');
end
fprintf(ID, '};\n');
[r,c]=size(IM);
fprintf(ID, '\n\n//Item Memory\nuint32_t iM[%d][%d] = {\n', r, c);
 
for i = 1 : r
    fprintf(ID, '{');
    for j = 1 : c
        if j == c
            fprintf(ID, '%d',(IM(i, j)));
        else 
            fprintf(ID, '%d, ',(IM(i, j)));  
        end
    end
    
    fprintf(ID, '},\n');
end
fprintf(ID, '};\n');

[r,c]=size(CIM);
fprintf(ID, '\n\n//Continuous Item Memory\nuint32_t ciM[%d][%d] = {\n', r, c);
 
for i = 1 : r
    fprintf(ID, '{');
    for j = 1 : c
        if j == c
            fprintf(ID, '%d',(CIM(i, j)));
        else 
            fprintf(ID, '%d, ',(CIM(i, j)));  
        end
    end
    
    fprintf(ID, '},\n');
end
fprintf(ID, '};\n');

[r,c]=size(AM);
fprintf(ID, '\n\n//Associative Memory\nuint32_t aM_32[%d][%d] = {\n', r, c);
 
for i = 1 : r
    fprintf(ID, '{');
    for j = 1 : c
        if j == c
            fprintf(ID, '%d',(AM(i, j)));
        else 
            fprintf(ID, '%d, ',(AM(i, j)));  
        end
    end
    
    fprintf(ID, '},\n');
end

fprintf(ID, '};');
fprintf(ID, '\n\n#endif');
fclose(ID);