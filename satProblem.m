classdef satProblem < handle
    %SATPROBLEM Summary of this class goes here
    %   Detailed explanation goes here

    properties
        decisionVariables
        nDecisionVariables
        clauses
    end

    methods
        function v = var(obj, varName, varargin)
            s = obj.decisionVariables(varName);
            s = s{1};
            v = s(varargin{:});
        end
        function obj = satProblem()
            %SATPROBLEM Construct an instance of this class
            %   Detailed explanation goes here
            %obj.Property1 = inputArg1 + inputArg2;
            obj.decisionVariables = dictionary();
            obj.nDecisionVariables = 0;
        end

        function obj = addDecisionVariable(obj, varName, mask)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            first = obj.nDecisionVariables + 1; %starting index
            last = obj.nDecisionVariables + sum(mask, 'all'); %ending index
            obj.nDecisionVariables = last;

            N = zeros(size(mask));

            N(logical(mask(:))) = first:1:last;

            obj.decisionVariables(varName) = {N};
        end
        function obj = addClause(obj, vec)
            obj.clauses{end+1} = vec;
        end

        function s = isSat(obj, solver)

            if strcmp(solver, 'minisat')
                problemName = obj.writeProblem(inputname(1));
                command = ['minisat ', [problemName, '.in '], [problemName, '.out']];

                path1 = getenv('PATH');
                if isempty(strfind(path1, ['/opt/local/bin/', pathsep]))
                    path1 = [path1 ':/opt/local/bin/'];
                    setenv('PATH', path1);
                end

                if exist([problemName, '.out'], 'file')==2
                    delete([problemName, '.out']);
                end


                [status,cmdout] = system(command);
                %display(status)
                %display(cmdout)

            elseif strcmp(solver, 'glucose-syrup')
                problemName = obj.writeProblem(inputname(1));
                command = ['glucose-syrup ', [problemName, '.in '], [problemName, '.out']];

                path1 = getenv('PATH');
                if isempty(strfind(path1, ['/Applications/glucose-main/parallel/', pathsep]))
                    path1 = [path1 ':/Applications/glucose-main/parallel/'];
                    setenv('PATH', path1);
                end

                if exist([problemName, '.out'], 'file')==2
                    delete([problemName, '.out']);
                end


                [status,cmdout] = system(command);
                %display(status)
                %display(cmdout)

            elseif strcmp(solver, 'glucose')
                problemName = obj.writeProblem(inputname(1));
                command = ['glucose -verb=1 ', [problemName, '.in '], [problemName, '.out']];

                path1 = getenv('PATH');
                if isempty(strfind(path1, ['/Applications/glucose-main/simp/', pathsep]))
                    path1 = [path1 ':/Applications/glucose-main/simp/'];
                    setenv('PATH', path1);
                end

                if exist([problemName, '.out'], 'file')==2
                    delete([problemName, '.out']);
                end
                fid = fopen([problemName, '.out'], 'wt');
                fclose(fid);


                [status,cmdout] = system(command);
                %display(status)
                %display(cmdout)

            end

            S = readlines([problemName, '.out']);

            if strcmp(S(1), 'SAT')
                %satisfiable
                s.sat = 1;

                %get assignment
                s.assignment = str2num(S(2));
                s.assignment(end) = [];

                %minisat doesn't care how many decision variables there
                %are, just how many
                s.values = dictionary();
                k = keys(obj.decisionVariables, 'cell');
                for i = 1:length(k)
                    I = obj.decisionVariables(k{i});
                    I = I{1};
                    M = -ones(size(I));

                    V = ( abs(s.assignment) >= min(I(I(:) > 0)) ) & ( abs(s.assignment) <= max(I(I(:)>0)) );

                    M(I(:) > 0) = s.assignment(V) > 0;

                    s.values(k{i}) = {M};
                end
            elseif strcmp(S(1), 'UNSAT')
                %unsatisfiable
                s.sat = 0;
            else
                %oops
                display(S(1))
                error('unknown exit condition of solver')
            end


        end
    end
    methods (Access = private)
        function problemName = writeProblem(obj, problemName)
            %creates a DIMACS file for the SAT problem
            %problemName = inputname(1);

            h = fopen([problemName, '.in'], 'w');
            fprintf(h, ['p cnf ', num2str(obj.nDecisionVariables), ' ', num2str(length(obj.clauses)), '\n']);

            for i = 1:length(obj.clauses)
                fprintf(h, '%d ', obj.clauses{i});
                fprintf(h, '0');
                if i < length(obj.clauses)
                    fprintf(h, '\n');
                end
            end

            fclose(h);
        end
    end
end

