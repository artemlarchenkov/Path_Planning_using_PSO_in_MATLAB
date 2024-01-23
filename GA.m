clc;
clear;
close all;

%% Определение задачи

model = CreateModel();  % Создание структуры, представляющей модель

model.n = 3;  % Количество точек управления

CostFunction = @(x) MyCost(x, model);  % Функция стоимости

nVar = model.n;       % Количество переменных решения

VarSize = [1, nVar];   % Размер матрицы переменных решения

VarMin.x = model.xmin;           % Нижняя граница переменных
VarMax.x = model.xmax;           % Верхняя граница переменных
VarMin.y = model.ymin;           % Нижняя граница переменных
VarMax.y = model.ymax;           % Верхняя граница переменных

%% Параметры генетического алгоритма

MaxIt = 500;          % Максимальное количество итераций

nPop = 100;           % Размер популяции (количество индивидуумов)

pc = 0.9;             % Вероятность кроссинговера
nc = 2 * round(pc * nPop / 2);  % Количество индивидуумов для кроссинговера (четное)

pm = 1;             % Вероятность мутации
nm = round(pm * nPop);  % Количество индивидуумов для мутации

%% Инициализация популяции

empty_particle.Position = struct('x', [], 'y', []);
empty_particle.Cost = [];
empty_particle.Sol = [];

pop = repmat(empty_particle, nPop, 1);

for i = 1:nPop
    pop(i).Position = CreateRandomSolution(model);
    [pop(i).Cost, pop(i).Sol] = CostFunction(pop(i).Position);
end

% Массив для хранения лучших значений стоимости на каждой итерации
BestCostGA = zeros(MaxIt, 1);

UpdateFrequency = 5;  % Частота обновления графика

%% Основной цикл генетического алгоритма

for it = 1:MaxIt
    
    % Оценка приспособленности
    costs = [pop.Cost];
    [costs, sortedIndices] = sort(costs);
    pop = pop(sortedIndices);
    
    % Сохранение лучшего решения
    BestSolGA = pop(1).Sol;
    BestCostGA(it) = pop(1).Cost;
    
    % Кроссинговер
    for k = 1:nc/2
        i1 = 2 * k - 1;
        i2 = 2 * k;
        [pop(i1).Position, pop(i2).Position] = Crossover(pop(i1).Position, pop(i2).Position);
        pop(i1).Position.x = max(VarMin.x, min(VarMax.x, pop(i1).Position.x));
        pop(i1).Position.y = max(VarMin.y, min(VarMax.y, pop(i1).Position.y));
        
        pop(i2).Position.x = max(VarMin.x, min(VarMax.x, pop(i2).Position.x));
        pop(i2).Position.y = max(VarMin.y, min(VarMax.y, pop(i2).Position.y));

    end
    
    % Мутация
    for k = 1:nm
        i = nPop - k + 1;
        pop(i).Position = Mutate(pop(i).Position, VarMin, VarMax);
        [pop(i).Cost, pop(i).Sol] = CostFunction(pop(i).Position);
    end
    
    % Обновление лучшей найденной стоимости
    BestCostGA(it) = pop(1).Cost;

    % Вывод информации о текущей итерации
    if  BestSolGA.IsFeasible
        Flag = ' *';
    else
        Flag = [', Нарушение = ' num2str( BestSolGA.Violation)];
    end
    disp(['Итерация ' num2str(it) ': Лучшая стоимость = ' num2str(BestCostGA(it)) Flag]);
    
    % Построение решения
    figure(2);
    PlotSolution( BestSolGA, model);
    axis([model.xmin model.xmax model.ymin model.ymax]);  % Установка пределов осей
    pause(0.01);

    % Построение решения с указанной частотой
    if mod(it, UpdateFrequency) == 0
        figure(2);
        PlotSolution(BestSolGA, model);
        pause(0.01);
    end
end

%% Результаты

figure;
plot(BestCostGA, 'LineWidth', 2);
xlabel('Итерация');
ylabel('Лучшая стоимость');
grid on;

function [child1, child2] = Crossover(parent1, parent2)
    nVar = numel(parent1.x);
    
    % Randomly select a crossover point
    crossoverPoint = randi([1, nVar - 1]);
    
    % Perform crossover
    child1.x = [parent1.x(1:crossoverPoint), parent2.x(crossoverPoint + 1:end)];
    child1.y = [parent1.y(1:crossoverPoint), parent2.y(crossoverPoint + 1:end)];
    
    child2.x = [parent2.x(1:crossoverPoint), parent1.x(crossoverPoint + 1:end)];
    child2.y = [parent2.y(1:crossoverPoint), parent1.y(crossoverPoint + 1:end)];
end

function mutatedPosition = Mutate(position, VarMin, VarMax)
    nVar = numel(position.x);
    
    % Randomly select a mutation point
    mutationPoint = randi([1, nVar]);
    
    % Generate a random value for mutation
    mutationValue = randn(size(position.x(mutationPoint)));
    
    % Apply mutation to the selected point
    mutatedPosition.x = position.x;
    mutatedPosition.x(mutationPoint) = position.x(mutationPoint) + mutationValue;
    
    mutationValue = randn(size(position.y(mutationPoint)));
    mutatedPosition.y = position.y;
    mutatedPosition.y(mutationPoint) = position.y(mutationPoint) + mutationValue;
    
    % Ensure mutated position is within bounds
    mutatedPosition.x = max(VarMin.x, min(VarMax.x, mutatedPosition.x));
    mutatedPosition.y = max(VarMin.y, min(VarMax.y, mutatedPosition.y));
end

