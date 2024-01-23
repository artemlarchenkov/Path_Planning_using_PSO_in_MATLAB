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

%% Параметры PSO

MaxIt = 500;          % Максимальное количество итераций

nPop = 150;           % Размер популяции (количество частиц)

w = 1;                % Инерционный вес
wdamp = 0.98;         % Скорость затухания инерционного веса
c1 = 1.5;             % Коэффициент личного обучения
c2 = 1.5;             % Коэффициент глобального обучения

alpha = 0.1;
VelMax.x = alpha * (VarMax.x - VarMin.x);    % Максимальная скорость
VelMin.x = -VelMax.x;                       % Минимальная скорость
VelMax.y = alpha * (VarMax.y - VarMin.y);    % Максимальная скорость
VelMin.y = -VelMax.y;                       % Минимальная скорость

%% Инициализация

% Создание пустой структуры для частиц
empty_particle.Position = struct('x', [], 'y', []);
empty_particle.Velocity = struct('x', [], 'y', []);
empty_particle.Cost = [];
empty_particle.Sol = [];
empty_particle.Best.Position = struct('x', [], 'y', []);
empty_particle.Best.Cost = [];
empty_particle.Best.Sol = [];

% Инициализация глобального лучшего решения
GlobalBest.Cost = inf;

% Создание матрицы частиц
particle = repmat(empty_particle, nPop, 1);

% Цикл инициализации
for i = 1:nPop
    
    % Инициализация положения частицы
    if i > 1
        particle(i).Position = CreateRandomSolution(model);
    else
        % Прямая линия от начала к концу
        xx = linspace(model.xs, model.xt, model.n + 2);
        yy = linspace(model.ys, model.yt, model.n + 2);
        particle(i).Position.x = xx(2:end-1);
        particle(i).Position.y = yy(2:end-1);
    end
    
    % Инициализация скорости
    particle(i).Velocity.x = zeros(VarSize);
    particle(i).Velocity.y = zeros(VarSize);
    
    % Оценка
    [particle(i).Cost, particle(i).Sol] = CostFunction(particle(i).Position);
    
    % Обновление лучшего решения частицы
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    particle(i).Best.Sol = particle(i).Sol;
    
    % Обновление глобального лучшего решения
    if particle(i).Best.Cost < GlobalBest.Cost
        GlobalBest = particle(i).Best;
    end
    
end

% Массив для хранения лучших значений стоимости на каждой итерации
BestCost = zeros(MaxIt, 1);

%% Основной цикл PSO

for it = 1:MaxIt
    
    for i = 1:nPop
        
        % x-часть
        
        % Обновление скорости
        particle(i).Velocity.x = w * particle(i).Velocity.x ...
            + c1 * rand(VarSize) .* (particle(i).Best.Position.x - particle(i).Position.x) ...
            + c2 * rand(VarSize) .* (GlobalBest.Position.x - particle(i).Position.x);
        
        % Ограничение скорости
        particle(i).Velocity.x = max(particle(i).Velocity.x, VelMin.x);
        particle(i).Velocity.x = min(particle(i).Velocity.x, VelMax.x);
        
        % Обновление положения
        particle(i).Position.x = particle(i).Position.x + particle(i).Velocity.x;
        
        % Зеркалирование скорости
        OutOfTheRange = (particle(i).Position.x < VarMin.x | particle(i).Position.x > VarMax.x);
        particle(i).Velocity.x(OutOfTheRange) = -particle(i).Velocity.x(OutOfTheRange);
        
        % Ограничение положения
        particle(i).Position.x = max(particle(i).Position.x, VarMin.x);
        particle(i).Position.x = min(particle(i).Position.x, VarMax.x);
        
        % y-часть
        
        % Обновление скорости
        particle(i).Velocity.y = w * particle(i).Velocity.y ...
            + c1 * rand(VarSize) .* (particle(i).Best.Position.y - particle(i).Position.y) ...
            + c2 * rand(VarSize) .* (GlobalBest.Position.y - particle(i).Position.y);
        
        % Ограничение скорости
        particle(i).Velocity.y = max(particle(i).Velocity.y, VelMin.y);
        particle(i).Velocity.y = min(particle(i).Velocity.y, VelMax.y);
        
        % Обновление положения
        particle(i).Position.y = particle(i).Position.y + particle(i).Velocity.y;
        
        % Зеркалирование скорости
        OutOfTheRange = (particle(i).Position.y < VarMin.y | particle(i).Position.y > VarMax.y);
        particle(i).Velocity.y(OutOfTheRange) = -particle(i).Velocity.y(OutOfTheRange);
        
        % Ограничение положения
        particle(i).Position.y = max(particle(i).Position.y, VarMin.y);
        particle(i).Position.y = min(particle(i).Position.y, VarMax.y);
        
        % Оценка
        [particle(i).Cost, particle(i).Sol] = CostFunction(particle(i).Position);
        
        % Обновление личного лучшего решения
        if particle(i).Cost < particle(i).Best.Cost
            
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;
            particle(i).Best.Sol = particle(i).Sol;
            
            % Обновление глобального лучшего решения
            if particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest = particle(i).Best;
            end
            
        end
        
    end
    
    % Обновление лучшей найденной стоимости
    BestCost(it) = GlobalBest.Cost;
    
    % Затухание инерционного веса
    w = w * wdamp;

    % Вывод информации о текущей итерации
    if GlobalBest.Sol.IsFeasible
        Flag = ' *';
    else
        Flag = [', Нарушение = ' num2str(GlobalBest.Sol.Violation)];
    end
    disp(['Итерация ' num2str(it) ': Лучшая стоимость = ' num2str(BestCost(it)) Flag]);
    
    % Построение решения
    figure(1);
    PlotSolution(GlobalBest.Sol, model);
    pause(0.01);    
end

%% Результаты

figure;
plot(BestCost, 'LineWidth', 2);
xlabel('Итерация');
ylabel('Лучшая стоимость');
grid on;

% Генетический алгоритм параметры
GA.MaxIt = 100;          % Максимальное количество итераций ГА
GA.nPop = 100;           % Размер популяции ГА
GA.pc = 0.8;             % Вероятность кроссовера
GA.gamma = 0.1;          % Вероятность мутации
GA.mu = 0.1;             % Вероятность изменения гена в мутации

% Инициализация популяции
pop = CreateParticles(GA.nPop, VarSize, VarMin, VarMax, model, CostFunction);

% Главный цикл ГА
BestCostGA = zeros(GA.MaxIt, 1);

[BestSolGA, BestCostGA] = GeneticAlgorithm(model, CostFunction, VarSize, VarMin, VarMax, GA);


function [BestSolGA, BestCostGA] = GeneticAlgorithm(model, CostFunction, VarSize, VarMin, VarMax, GA)

    % Инициализация популяции
    pop = CreateParticles(GA.nPop, VarSize, VarMin, VarMax, model, CostFunction);

    % Главный цикл ГА
    BestCostGA = zeros(GA.MaxIt, 1);

    for it = 1:GA.MaxIt

        % Вычисление значений стоимости для каждой особи
        costs = zeros(GA.nPop, 1);
        for i = 1:GA.nPop
            [costs(i), ~] = CostFunction(pop(i).Position);
        end

        TournamentSize = 3; % Выберите желаемый размер турнира
        parents = TournamentSelection(pop, [pop.Cost], TournamentSize);

        % Кроссовер
        offspring = Crossover(parents, GA.pc);

        % Мутация
        offspring = Mutate(offspring, GA.mu, GA.gamma, VarMin, VarMax);

        % Вычисление значений стоимости для потомства
        for i = 1:numel(offspring)
            [offspring(i).Cost, offspring(i).Sol] = CostFunction(offspring(i).Position);
        end

        % Объединение текущей популяции и потомства
        pop = [pop; offspring];

        % Сортировка популяции по стоимости
        [~, sortOrder] = sort([pop.Cost]);
        pop = pop(sortOrder);

        % Отбор лучших особей
        pop = pop(1:GA.nPop);

        % Обновление лучшей стоимости
        BestCostGA(it) = pop(1).Cost;

        % Вывод информации о текущей итерации
        disp(['GA Iteration ' num2str(it) ': Best Cost = ' num2str(BestCostGA(it))]);

    end

    % Постобработка результатов ГА
    BestSolGA = pop(1).Sol;

end


function particles = CreateParticles(nPop, VarSize, VarMin, VarMax, model, CostFunction)
    particles = struct('Position', [], 'Velocity', [], 'Cost', [], 'Sol', [], 'Best', struct('Position', [], 'Cost', [], 'Sol', []));

    for i = 1:nPop
        if i > 1
            particles(i).Position = CreateRandomSolution(model);
        else
            % Straight line from source to destination
            xx = linspace(model.xs, model.xt, model.n+2);
            yy = linspace(model.ys, model.yt, model.n+2);
            particles(i).Position.x = xx(2:end-1);
            particles(i).Position.y = yy(2:end-1);
        end

        particles(i).Velocity.x = zeros(VarSize);
        particles(i).Velocity.y = zeros(VarSize);

        [particles(i).Cost, particles(i).Sol] = CostFunction(particles(i).Position);

        particles(i).Best.Position = particles(i).Position;
        particles(i).Best.Cost = particles(i).Cost;
        particles(i).Best.Sol = particles(i).Sol;
    end
end

function parents = TournamentSelection(pop, costs, TournamentSize)
    nPop = numel(pop);
    parents = zeros(2, nPop);

    for i = 1:nPop
        % Выбираем случайные участников турнира
        tournamentInds = randperm(nPop, TournamentSize);
        
        % Находим победителя турнира (минимальная стоимость)
        [~, winnerInd] = min(costs(tournamentInds));
        
        % Выбираем победителя в качестве родителя
        parents(1, i) = tournamentInds(winnerInd);
    end

    % Выбираем второго родителя из оставшихся участников
    parents(2, :) = randi(nPop, 1, nPop);
end

function sol = CreateRandomSolution(model)
    sol.x = model.xmin + rand(1, model.n) * (model.xmax - model.xmin);
    sol.y = model.ymin + rand(1, model.n) * (model.ymax - model.ymin);
end
function offspring = Crossover(parents, pc)
    nParents = size(parents, 1);
    nVar = numel(parents(1).Position.x);

    offspring = struct('Position', struct('x', zeros(1, nVar), 'y', zeros(1, nVar)), ...
                      'Velocity', struct('x', zeros(1, nVar), 'y', zeros(1, nVar)), ...
                      'Cost', [], 'Sol', [], 'Best', struct('Position', struct('x', [], 'y', []), 'Cost', [], 'Sol', []));

    for i = 1:2:nParents
        % Выбор родителей
        parent1 = parents(i);
        parent2 = parents(i + 1);

        % Вероятность кроссовера
        if rand() <= pc
            % Выбираем случайную точку разрыва
            crossoverPoint = randi([1, nVar - 1]);

            % Обновляем потомков
            offspring(i).Position.x(1:crossoverPoint) = parent1.Position.x(1:crossoverPoint);
            offspring(i).Position.x(crossoverPoint + 1:end) = parent2.Position.x(crossoverPoint + 1:end);

            offspring(i + 1).Position.x(1:crossoverPoint) = parent2.Position.x(1:crossoverPoint);
            offspring(i + 1).Position.x(crossoverPoint + 1:end) = parent1.Position.x(crossoverPoint + 1:end);

            offspring(i).Position.y(1:crossoverPoint) = parent1.Position.y(1:crossoverPoint);
            offspring(i).Position.y(crossoverPoint + 1:end) = parent2.Position.y(crossoverPoint + 1:end);

            offspring(i + 1).Position.y(1:crossoverPoint) = parent2.Position.y(1:crossoverPoint);
            offspring(i + 1).Position.y(crossoverPoint + 1:end) = parent1.Position.y(crossoverPoint + 1:end);
        else
            % Если кроссовер не произошел, просто копируем родителей в потомков
            offspring(i) = parent1;
            offspring(i + 1) = parent2;
        end
    end
end


function pop = Mutate(pop, mu, gamma, VarMin, VarMax)
    nPop = numel(pop);
    nVar = numel(pop(1).Position.x);

    for i = 1:nPop
        % Мутация по каждой переменной
        for j = 1:nVar
            if rand() <= mu
                % Аддитивная мутация
                dx = gamma * (VarMax.x - VarMin.x);
                pop(i).Position.x(j) = pop(i).Position.x(j) + dx;

                dy = gamma * (VarMax.y - VarMin.y);
                pop(i).Position.y(j) = pop(i).Position.y(j) + dy;
            end
        end
    end
end