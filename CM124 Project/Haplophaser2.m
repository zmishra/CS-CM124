function Haplophaser2(input, output, partition_size)

%Parameters
SNP_Partition = partition_size;
DELTA = 10^-3;

%Load data
A = importdata(input, ' ');
phased_population = zeros(size(A, 1), 2*size(A, 2));

for start_ = 1:SNP_Partition:size(A,1)
    %Initialize and create partition
    end_ = min(start_+SNP_Partition-1, size(A,1));
    B = A(start_:end_,:);
    SNP_Count = size(B, 1);
    phased_partition = zeros(size(B, 1), 2*size(B, 2));
    
    %Find all haplotypes for the partition
    all_haplotypes = [];
    gen_haplotypes = cell(size(A, 2),1);
    gen_freq = cell(size(A, 2), 1);
    
    for i = 1:size(A, 2)
        hetero = [1 0];
        gen = zeros(SNP_Count,1);
        for j = 1:SNP_Count
            switch B(j, i)
                case 1
                    gen = [gen gen];
                    if size(hetero, 2) < size(gen,2)
                        hetero = [hetero ~hetero];
                    end
                    gen(j, :) = hetero;
                case 2
                    gen(j,:) = 1;
                case 0
                    gen(j,:) = 0;
            end
        end
        gen_haplotypes{i} = gen;
        gen_freq{i} = ones(1, ceil(size(gen, 2)/2))/ceil(size(gen, 2)/2);

        all_haplotypes = unique([all_haplotypes gen].', 'rows').';
    end
    
    %Create a map for easy look-up of haplotype index
    allhapcells = cell(1, size(all_haplotypes, 2));
    for i = 1:1:size(all_haplotypes, 2)
        allhapcells{i} = num2str(all_haplotypes(:,i).');
    end
    hapMap = containers.Map(allhapcells, 1:1:size(all_haplotypes,2));
    
    %Initialize haplotype frequencies
    P = zeros(size(all_haplotypes,2),1);
    for i = 1:size(gen_haplotypes,1)
        P_gen = zeros(size(all_haplotypes,2),1);
        curr_gen = gen_haplotypes{i};
        curr_gen_freq = gen_freq{i};
%         P_max = 0;
        
        if size(curr_gen, 2) ~= 1
            for j = 1:2:size(curr_gen, 2)
                hap1 = curr_gen(:,j);
                hap2 = curr_gen(:,j+1);
                index_hap1 = hapMap(num2str(hap1.'));
                index_hap2 = hapMap(num2str(hap2.'));
                P_gen(index_hap1) = curr_gen_freq((j+1)/2);
                P_gen(index_hap2) = P_gen(index_hap1);
            end
        else
            index_hap1 = hapMap(num2str(curr_gen(:,1).'));
            P_gen(index_hap1) = 2;
        end
        P = P + P_gen;
    end
    P = P/size(gen_haplotypes,1)/2;
%   P = 1/size(all_haplotypes,2)*ones(size(all_haplotypes,2),1);
    P_old = ones(size(all_haplotypes,2),1)*-1000;
    index_P_max = zeros(size(gen_haplotypes,1),2);
    
    iterations = 0;
    blacklist=zeros(size(gen_haplotypes,1),1);

    %EM
    while max(abs(P-P_old)) > DELTA && sum(blacklist) ~= size(blacklist, 1)
        iterations = iterations + 1;
        
        
        P_new = zeros(size(all_haplotypes,2),1);
        for i = 1:size(gen_haplotypes,1)
            P_gen = zeros(size(all_haplotypes,2),1);
            if blacklist(i) == 0
            curr_gen = gen_haplotypes{i};
            P_max = 0;
            
            if size(curr_gen, 2) ~= 1
                for j = 1:2:size(curr_gen, 2)
                    hap1 = curr_gen(:,j);
                    hap2 = curr_gen(:,j+1);
                    index_hap1 = hapMap(num2str(hap1.'));
                    index_hap2 = hapMap(num2str(hap2.'));
                    P_gen(index_hap1) = P(index_hap1)*P(index_hap2);
                    P_gen(index_hap2) = P_gen(index_hap1);
                    if P_gen(index_hap1) > P_max
                        P_max = P_gen(index_hap1);
                        index_P_max(i,1) = index_hap1;
                        index_P_max(i,2) = index_hap2;
                    end
                end
                P_gen = P_gen/sum(P_gen);
                P_gen(P_gen < 10^-2) = 0;
                P_gen = P_gen/sum(P_gen);
                if max(P_gen)>0.25
                    blacklist(i) = 1;
                end
            else
                index_hap1 = hapMap(num2str(curr_gen(:,1).'));
                P_gen(index_hap1) = 1;
                index_P_max(i,1) = index_hap1;
                index_P_max(i,2) = index_hap1;
                blacklist(i) = 1;
            end
            else
                P_gen(index_P_max(i,1)) = P_gen(index_P_max(i,1)) + 0.5;
                P_gen(index_P_max(i,2)) = P_gen(index_P_max(i,2)) + 0.5;
            end
            P_new = P_new + P_gen;
        end
        P_old = P;
        
        
        P = P_old + 1*(P_new/size(gen_haplotypes,1)-P_old);
    end
    
    %Construct phased partition
    for i = 1:size(gen_haplotypes,1)
        phased_partition(:,2*i-1) = all_haplotypes(:,index_P_max(i,1));
        phased_partition(:,2*i) = all_haplotypes(:,index_P_max(i,2));
    end
    
    %Add phased partition to phased population
    phased_population(start_:end_,:) = phased_partition;
    
end

dlmwrite(output, phased_population, 'delimiter', ' ');



