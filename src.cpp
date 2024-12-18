#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
#include <algorithm>
class Line { // линия представляет собой две точки её концов
public:
	double x1;
	double y1;
	double x2;
	double y2;
	static bool checkDot(double x, double x1, double x2) { // проверка лежит ли x между x1 и x2
		if (x < std::max(x1, x2) && x > std::min(x1, x2)) {
			return true;
		}
		else {
			return false;
		}
	}
	Line() {}
	Line(double _x1, double _y1, double _x2, double _y2) : x1(_x1), x2(_x2), y1(_y1), y2(_y2) {};
	bool intersect(Line line) { //проверка пересекаются ли данная линия с линией line
		int mod1 = x2 - x1;
		int mod2 = line.x2 - line.x1;
		if ((x1 == line.x1 && y1 == line.y1) || (x1 == line.x2 && y1 == line.y2) || //проверка если рёбра имеют общую вершину
			(x2 == line.x1 && y2 == line.y1) || (x2 == line.x2 && y2 == line.y2)) {
			if (mod1 == 0) {
				if (mod2 == 0) {
					if (x1 == line.x1 && y1 == line.y1 ) {
						return checkDot(x2, line.x1, line.x2) || checkDot(line.x2, x1, x2);
					}
					else if (x1 == line.x2 && y1 == line.y2) {
						return checkDot(x2, line.x1, line.x2) || checkDot(line.x1, x1, x2);
					}
					else if (x2 == line.x1 && y2 == line.y1) {
						checkDot(x1, line.x1, line.x2) || checkDot(line.x2, x1, x2);
					}
					else {
						checkDot(x1, line.x1, line.x2) || checkDot(line.x1, x1, x2);
					}
				}
				else {
					return false;
				}
			}
			else {
				if (mod2 == 0) {
					return false;
				}
				else {
					double k1 = ((y2 - y1) * 1.0) / mod1;
					double k2 = ((line.y2 - line.y1) * 1.0) / mod2;
					if (k1 != k2) {
						return false;
					}
					else {
						if (x1 == line.x1 && y1 == line.y1) {
							return checkDot(x2, line.x1, line.x2) || checkDot(line.x2, x1, x2);
						}
						else if (x1 == line.x2 && y1 == line.y2) {
							return checkDot(x2, line.x1, line.x2) || checkDot(line.x1, x1, x2);
						}
						else if (x2 == line.x1 && y2 == line.y1) {
							checkDot(x1, line.x1, line.x2) || checkDot(line.x2, x1, x2);
						}
						else {
							checkDot(x1, line.x1, line.x2) || checkDot(line.x1, x1, x2);
						}
					}
				}
			}
		}
		if (mod1 == 0) {
			if (mod2 == 0) {
				if (x1 == line.x1) {
					return checkDot(line.y1, y1, y2) || checkDot(line.y2, y1, y2);
				}
				else {
					return false;
				}
			}
			else {
				double k2 = ((line.y2 - line.y1) * 1.0) / mod2;
				double b2 = line.y1 - k2 * line.x1;
				double y = k2 * x1 + b2;
				return checkDot(y, y1, y2) && checkDot(x1, line.x1, line.x2);
			}
		}
		else {
			double k1 = ((y2 - y1) * 1.0) / mod1;
			double b1 = y1 - k1 * x1;
			if (mod2 == 0) {
				double y = k1 * line.x1 + b1;
				return checkDot(y, line.y1, line.y2) && checkDot(line.x1, x1, x2);
			}
			else {
				double k2 = ((line.y2 - line.y1) * 1.0) / mod2;
				double b2 = line.y1 - k2 * line.x1;
				if (k2 != k1) {
					double x = (b2 - b1) / (k1 - k2);
					return checkDot(x, x1, x2) && checkDot(x, line.x1, line.x2);
				}
				else {
					if (b1 == b2) {
						return checkDot(line.x1, x1, x2) || checkDot(line.x2, x1, x2);
					}
					else {
						return false;
					}
				}
			}
		}
	}
};
struct Dot { // вершина, представляющая собой точку на плоскости
	Dot() {}
	Dot(double _x, double _y) : x(_x), y(_y) {}
	double x;
	double y;
};
struct Edge { //ребро, представляющее собой номера вершин 
	Edge() {}
	Edge(int _dot1, int _dot2) : dot1(_dot1), dot2(_dot2) {}
	int dot1;
	int dot2;
};
class Graph {
public:
	int verticesNum;
	int edgeNum;
	std::vector<Dot> vertices;
	std::vector<Edge> edges;
	std::vector<std::vector<int>> crossings;
	Graph(int vertNum, int edNum, std::vector<Dot> _vertices, std::vector<Edge> _edges) :verticesNum(vertNum), edgeNum(edNum),
		vertices(_vertices), edges(_edges) {}
	void getCrossing() {
		crossings = std::vector<std::vector<int>>(edgeNum);
		Line* lines = new Line[edgeNum];
		for (int i = 0; i < edgeNum; i++) {
			lines[i] = Line(vertices[edges[i].dot1].x, vertices[edges[i].dot1].y,
				vertices[edges[i].dot2].x, vertices[edges[i].dot2].y);
		}
		for (int i = 0; i < edgeNum; i++) { //проверка пересекаются ли отрезки
			for (int j = i + 1; j < edgeNum; j++) {
				if (lines[i].intersect(lines[j])) {
					crossings[i].push_back(j);
					crossings[j].push_back(i);
				}
			}
		}
		delete[] lines;
	}
	Graph getSubgraph(std::vector<int> subEdge) { //нахождения подграфа по списку рёбер
		std::vector<Edge> subEdges = std::vector<Edge>(subEdge.size());
		for (int i = 0; i < subEdge.size(); i++) {
			subEdges[i] = edges[subEdge[i]];
		}
		return Graph(verticesNum, subEdges.size(), vertices, subEdges);
	}
	void saveToTextFile(std::string name) { //сохранение графа в текстовый файл
		std::ofstream fout(name);
		fout << verticesNum << '\n';
		fout << edgeNum << '\n';
		for (int i = 0; i < verticesNum; i++) {
			fout << vertices[i].x << ' ' << vertices[i].y << '\n';
		}
		for (int i = 0; i < edgeNum; i++) {
			fout << edges[i].dot1 << ' ' << edges[i].dot2 << '\n';
		}
	}
	static Graph readFromFile(std::string name) { //чтение графа из текстового файла
		std::ifstream fin(name);
		int n, m;
		fin >> n >> m;
		std::vector<Dot> vert = std::vector<Dot>(n);
		std::vector<Edge> edge = std::vector<Edge>(m);
		for (int i = 0; i < n; i++) {
			int x, y;
			fin >> x >> y;
			vert[i] = Dot(x, y);
		}
		for (int i = 0; i < m; i++) {
			int dot1, dot2;
			fin >> dot1 >> dot2;
			edge[i] = Edge(dot1, dot2);
		}
		return Graph(n, m, vert, edge);
	}
};
void drawGraph(std::string fileName, std::string name, std::string color) { //рисование графа
	std::ostringstream command; //команда для вызова питоновского крипта
	command << "python3 graphVis.py ";
	command << fileName<<' ';
	command << color << ' ';
	command << name << ' ';
	system(command.str().c_str());
}
Graph createRandomGraph(int n, int m) {
	std::vector<Dot> vertices = std::vector<Dot>(n); //множество вершин
	int j = 0;
	while (j < n) { //добавление вершин
		bool flag = true;
		int x1 = rand() % 101;
		int y1 = rand() % 101;
		for (int i = 0; i < j; i++) {
			if (vertices[i].x == x1 && vertices[i].y == y1) { //проверка, что точка с данными координатами не существует
				flag = false;
				break;
			}
		}
		if (flag) {
			vertices[j] = Dot(x1, y1);
			j++;
		}
	}
	//генерация рёбер
	if (m < n * (n - 1) / 4) { 
		std::vector<Edge> edge = std::vector<Edge>(m);
		j = 0;
		while (j < m) {
			bool flag = true;
			int dot1 = rand() % n;
			int dot2 = rand() % n;
			if (dot1 == dot2) {
				continue;
			}
			for (int i = 0; i < j; i++) { // проверка что два ребра, не объединяли одни и те же вершины
				if (edge[i].dot1 == dot1 && edge[i].dot2 == dot2 || edge[i].dot2 == dot1 && edge[i].dot1 == dot2) {
					flag = false;
					break;
				}
			}
			if (flag) {
				edge[j] = Edge(dot1, dot2);
				j++;
			}
		}
		return Graph(n, m, vertices, edge);
	}
	else { //при большом количестве рёбер, рёбра выбираются случайным образом для дополнения графа
		int stop = n * (n - 1) / 2 - m;
		j = 0;
		std::vector<std::vector<bool>> matrixEdge = std::vector<std::vector<bool>>(n);
		for (int i = 0; i < n; i++) {
			matrixEdge[i] = std::vector<bool>(n);
			for (int k = 0; k < n; k++) {
				matrixEdge[i][k] = true;
			}
		}
		while (j < stop) { // проверка что два ребра, не объединяли одни и те же вершины
			bool flag = true;
			int dot1 = rand() % n;
			int dot2 = rand() % n;
			if (dot1 == dot2) {
				continue;
			}
			if (matrixEdge[dot1][dot2]) {
				matrixEdge[dot1][dot2] = false;
				matrixEdge[dot2][dot1] = false;
				j++;
			}
		}
		int it = 0;
		std::vector<Edge> edge = std::vector<Edge>(m);
		for (int i = 0; i < matrixEdge.size(); i++) {
			for (int k = i+1; k < matrixEdge.size(); k++) {
				if (matrixEdge[i][k]) {
					edge[it] = Edge(i, k);
					it++;
				}
			}
		}
		return Graph(n, m, vertices, edge);
	}
}
Graph createRandomCircleGraph(int vertNum, int edgeNum, int dist) {
	std::vector<Dot> vertices = std::vector<Dot>(vertNum); //множество вершин
	std::vector<int> distance = std::vector<int>(vertNum); //расстояние от соседа
	int resultDistance = 0; //общая сумма расстояний
	for (int i = 0; i < vertNum; i++) {
		distance[i] = 2*(rand() % dist+dist/2); //генерация расстояния
		resultDistance += distance[i];
	}
	double circleAngle = 0; //угол который соответствует вершине в круге
	double r = resultDistance / 2.0; // радиус круга
	for (int i = 0; i < vertNum; i++) {
		vertices[i] = Dot(r * cos(circleAngle)+resultDistance, r * sin(circleAngle) + 2*r + resultDistance);
		circleAngle += 2 * M_PI*distance[i] / (resultDistance);
	}
	if (edgeNum > vertNum * (vertNum - 1) / 4) { //генерация рёбер аналогичная createRandomGraph
		std::vector<Edge> edges = std::vector<Edge>(edgeNum);
		int currentEdgeNum = 0;
		while (currentEdgeNum < edgeNum) {
			bool flag = true;
			int dot1 = rand() % vertNum;
			int dot2 = rand() % vertNum;
			if (dot1 == dot2) {
				continue;
			}
			for (int i = 0; i < currentEdgeNum; i++) {
				if (edges[i].dot1 == dot1 && edges[i].dot2 == dot2 || edges[i].dot2 == dot1 && edges[i].dot1 == dot2) {
					flag = false;
					break;
				}
			}
			if (flag) {
				edges[currentEdgeNum] = Edge(dot1, dot2);
				currentEdgeNum++;
			}
		}
		return Graph(vertNum, edgeNum, vertices, edges);
	}
	else {
		int stopEdgeNum = vertNum * (vertNum - 1) / 2 - edgeNum;
		int currentEdgeNum = 0;
		std::vector<std::vector<bool>> matrixEdge = std::vector<std::vector<bool>>(vertNum);
		for (int i = 0; i < vertNum; i++) {
			matrixEdge[i] = std::vector<bool>(vertNum);
			for (int k = 0; k < vertNum; k++) {
				matrixEdge[i][k] = true;
			}
		}
		while (currentEdgeNum < stopEdgeNum) {
			bool flag = true;
			int dot1 = rand() % vertNum;
			int dot2 = rand() % vertNum;
			if (dot1 == dot2) {
				continue;
			}
			if (matrixEdge[dot1][dot2]) {
				matrixEdge[dot1][dot2] = false;
				matrixEdge[dot2][dot1] = false;
				currentEdgeNum++;
			}
		}
		int it = 0;
		std::vector<Edge> edge = std::vector<Edge>(edgeNum);
		for (int i = 0; i < matrixEdge.size(); i++) {
			for (int k = i + 1; k < matrixEdge.size(); k++) {
				if (matrixEdge[i][k]) {
					edge[it] = Edge(i, k);
					it++;
				}
			}
		}
		return Graph(vertNum, edgeNum, vertices, edge);
	}
}
Graph createRandomMatroidGraph(int sqRow, int sqCol, int chance1, int chance2) {
	std::vector<Dot> vertices = std::vector<Dot>();
	std::vector<Edge> edges = std::vector<Edge>();
	for (int i = 0; i < sqRow; i++) {
		for (int j = 0; j < sqCol; j++) {
			vertices.push_back(Dot(0 + i * 20, 0 + j * 20));
			vertices.push_back(Dot(5 + i * 20, 0 + j * 20));
			vertices.push_back(Dot(10 + i * 20, 0 + j * 20));
			vertices.push_back(Dot(10 + i * 20, 5 + j * 20));
			vertices.push_back(Dot(10 + i * 20, 10 + j * 20));
			vertices.push_back(Dot(5 + i * 20, 10 + j * 20));
			vertices.push_back(Dot(0 + i * 20, 10 + j * 20));
			vertices.push_back(Dot(0 + i * 20, 5 + j * 20));
			for (int u = 0; u < 8; u++) {
				if (rand() % chance2 <= chance1) {
					edges.push_back(Edge(u + 8 * j + sqCol * i * 8, (u + 1) % 8 + 8 * j + sqCol * i * 8));
				}
			}
			for (int u = 0; u < 4; u++) {
				if (rand() % chance2 <= chance1) {
					edges.push_back(Edge(u + 8 * j + sqCol * i * 8, (u + 4) % 8 + 8 * j + sqCol * i * 8));
				}
			}
			if (i > 0) {
				if (rand() % chance2 <= chance1) {
					edges.push_back(Edge(2 + 8 * j + sqCol * (i-1) * 8,  8 * j + sqCol * i * 8));
				}
			}
			if (j > 0) {
				if (rand() % chance2 <= chance1) {
					edges.push_back(Edge(4 + 8 * (j-1) + sqCol * i * 8, 2 + 8 * j + sqCol * i * 8));
				}
			}
		}
	}
	std::cout << vertices.size() << ' ' << edges.size() << '\n';
	return Graph(vertices.size(), edges.size(), vertices, edges);
}
class DSU {
private:
	std::vector<int> set; //множество элементов
	std::vector<std::vector<int>> memorySet; // список изменений в СНП
	int memoryNum = 0; //количество пересечений множество, состояние памяти СНП
public:
	DSU(int n) {
		set = std::vector<int>(n);
		memorySet = std::vector<std::vector<int>>(n);
		for (int i = 0; i < n; i++) {
			set[i] = -1;
			memorySet[i] = std::vector<int>(4);
		}
	}
	int getParrent(int dot) { //Найти самого глубокого предка 
		while (set[dot] >= 0) {
			dot = set[dot];
		}
		return dot;
	}
	bool connect(int dot1, int dot2) { //объединить элементов в одно множество, возвращает true когда элементы были объединены, false когда элементы до этого были в одном множестве
		int parent1 = getParrent(dot1);
		int parent2 = getParrent(dot2);
		if (parent1 == parent2) {
			return false;
		}
		else {
			dot1;
			dot2;
			memorySet[memoryNum][0] = parent1;
			memorySet[memoryNum][1] = set[parent1];
			memorySet[memoryNum][2] = parent2;
			memorySet[memoryNum][3] = set[parent2];
			memoryNum++;
			if (set[parent1] < set[parent2]) {
				set[parent1] += set[parent2];
				set[parent2] = parent1;
			}
			else {
				set[parent2] += set[parent1];
				set[parent1] = parent2;
			}
			return true;
		}
	}
	void restoreMemory() { //Возвращает к предыдущему состоянию СНП
		if (memoryNum == 0) {
			std::cout << "what";
		}
		int parent1 = memorySet[memoryNum - 1][0];
		int parent2 = memorySet[memoryNum - 1][2];
		set[parent1] = memorySet[memoryNum - 1][1];
		set[parent2] = memorySet[memoryNum - 1][3];
		memoryNum--;
	}
	bool checkEdge(int num1, int num2) { //Проверка имеют элементы одного предка
		return getParrent(num1) == getParrent(num2);
	}
};
int calculateComponents(Graph graph) {
	DSU graphCycleDSU=DSU(graph.verticesNum); //СНП для проверки циклов
	int forestEdgeNum = 0; // количество рёбер
	for (int i = 0; i < graph.edgeNum; i++) {
		if (graphCycleDSU.connect(graph.edges[i].dot1, graph.edges[i].dot2)) {
			forestEdgeNum++;
		}
	}
	return graph.verticesNum - forestEdgeNum;
}
Graph baseAlgorithm(Graph graph) {
	DSU graphCyclesDSU = DSU(graph.verticesNum); //СНП для исключения циклов
	int maxEdgeNum = graph.verticesNum - calculateComponents(graph); //максимальное число рёбер в непересекающемся подграфе
	int resultEdgeIter = 0;
	std::vector<int> resultEdges;
	for (int i = 0; i < graph.edgeNum; i++) { //заполнение графа непересекающимися рёбрами
		if (graph.crossings[i].size() == 0) {
			if (graphCyclesDSU.connect(graph.edges[i].dot1, graph.edges[i].dot2)) {
				resultEdges.push_back(i);
				resultEdgeIter++;
			}
		}
	}
	std::vector<int> remainEdges = std::vector<int>();
	for (int i = 0; i < graph.edgeNum; i++) {
		if (!graphCyclesDSU.checkEdge(graph.edges[i].dot1, graph.edges[i].dot2)) { //удаление рёбер образующих цикл с непересекающимися
			remainEdges.push_back(i);
		}
	}
	int edgeIter = 0;
	int edgeIterNum = 0;
	int maxIterVal = maxEdgeNum - resultEdgeIter;
	std::vector<int> currentEdges = std::vector<int>(maxIterVal); // набор рёбер на нынешней итерации
	std::vector<int> maxEdges = std::vector<int>(maxIterVal); //максимально возможно найденный набор рёбер
	std::vector<std::vector<bool>> currentIntersectEdges = std::vector<std::vector<bool>>(maxIterVal+1); // список рёбер, пересекающийся с выбранными рёбрами
	for (int i = 0; i < currentIntersectEdges.size(); i++) {
		currentIntersectEdges[i] = std::vector<bool>(graph.edgeNum);
	}
	int maxEdgeIterNum = 0;
	while (edgeIterNum<maxIterVal) {
		if (edgeIter == remainEdges.size()) {
			if (edgeIterNum == 0) {
				break; //окончание цикла, когда мы пытаемся убрать ребро из пустого графа
			}
			else {
				graphCyclesDSU.restoreMemory(); //возвращение к предыдущему состоянию
				edgeIter = currentEdges[edgeIterNum - 1]+1;
				edgeIterNum--;
				continue;
			}
		}
		if (!currentIntersectEdges[edgeIterNum][remainEdges[edgeIter]]) { //проверка на пересечения
			if (graphCyclesDSU.connect(graph.edges[remainEdges[edgeIter]].dot1, graph.edges[remainEdges[edgeIter]].dot2)) { // проверка на цикл
				currentEdges[edgeIterNum] = edgeIter; //добавления ребра
				edgeIterNum++;
				for (int i = 0; i < graph.edgeNum; i++) { // копирование списка пересечений
					currentIntersectEdges[edgeIterNum][i] = currentIntersectEdges[edgeIterNum-1][i];
				}
				for (int i = 0; i < graph.crossings[remainEdges[edgeIter]].size(); i++) { //пополнение списка пересечений
					currentIntersectEdges[edgeIterNum][graph.crossings[remainEdges[edgeIter]][i]] = true;
				}
				if (edgeIterNum > maxEdgeIterNum) { // замена максимального набора рёбер
					maxEdgeIterNum = edgeIterNum;
					for (int i = 0; i < maxEdgeIterNum; i++) {
						maxEdges[i] = currentEdges[i];
					}
				}
			}
		}
		edgeIter++;
	}
	for (int i = 0; i < maxEdgeIterNum; i++) { // добавление рёбер, полученных в результате итераций
		resultEdges.push_back(remainEdges[maxEdges[i]]);
	}
	return graph.getSubgraph(resultEdges);
}
Graph convexAlgorithm(Graph graph) {
	std::vector<std::vector<int>> sectorMaxEdgeNum = std::vector<std::vector<int>>(graph.verticesNum); //максимальнок количество рёбер в неперескающемся подграфе
	std::vector<std::vector<int>> sectorPointUnion = std::vector<std::vector<int>>(graph.verticesNum); //вершина, по которой произошло объединение в подграф
	std::vector<std::vector<bool>> edgeMatrix = std::vector<std::vector<bool>>(graph.verticesNum); // матрица смежности
	std::vector<std::vector<bool>> isSectorWayExist = std::vector<std::vector<bool>>(graph.verticesNum); //содержит ли сектор путь соединяющий его края
	std::vector<std::vector<bool>> isEdgeAdd = std::vector<std::vector<bool>>(graph.verticesNum); // было ли добавлено ребро при объдинении графов
	/*Идея для нумерации такова, предполагается, что граф выпуклый,берутся три точки, на них строится треугольник, и берётся точка внутри треугольника,
	в данном случае середина средней линии треугольника, строится отрезки из данной точки ко всем точкам многоугольника и высчитывается соотношение координат x и y 
	данных отрезков, отчего оно сортируется таким образом: всё что левее точки, меньше всего что, правее точки, точка находящаяся выше, другой точки больше,
	если они находятся по одну сторону от точки*/

	std::vector<std::pair<int,Dot>> linesFromDot(graph.verticesNum); 
	if (graph.verticesNum > 2) {
		Dot vert1 = graph.vertices[0];
		Dot vert2 = graph.vertices[0];
		Dot vert3 = graph.vertices[0];
		double x1 = (vert1.x + vert2.x) / 2;
		double y1 = (vert1.y + vert2.y) / 2;
		double x2 = (vert1.x + vert3.x) / 2;
		double y2 = (vert1.y + vert3.y) / 2;
		Dot vert = Dot((x1 + x2) / 2, (y1 + y2) / 2); //внутренния точка
		for (int i = 0; i < graph.verticesNum; i++) {
			double x = graph.vertices[i].x - vert.x;
			double y = graph.vertices[i].y - vert.y;
			double len = sqrt(x * x + y * y);
			linesFromDot[i] = std::pair<int,Dot>(i,Dot(x/len,y/len)); //нормированные отрезки
		}
	}
	else {
		return graph; // граф из одной или двух вершин всегда не имеет пересечений
	}
	std::sort(linesFromDot.begin(), linesFromDot.end(), [](std::pair<int, Dot> x, std::pair<int, Dot> y) { // сортировка отрезков
		if (x.second.x > 0) {
			if (y.second.x > 0) {
				return x.second.y > y.second.y;
			}
			else {
				return true;
			}
		}
		else if (y.second.x>0){
			return false;
		}
		else {
			return x.second.y < y.second.y;
		}
		});
	std::vector<int> enumerateVertices(graph.verticesNum); //извлечение отсортированных по отрезкам вершин
	for (int i = 0; i < graph.verticesNum; i++) {
		enumerateVertices[i] = linesFromDot[i].first;
	}

	for (int i = 0; i < graph.verticesNum;i++) { //заполнение всех матриц значениями, по умолчанию
		edgeMatrix[i] = std::vector<bool>(graph.verticesNum);
		sectorMaxEdgeNum[i] = std::vector<int>(graph.verticesNum,0);
		sectorPointUnion[i] = std::vector<int>(graph.verticesNum,-1);
		isSectorWayExist[i] = std::vector<bool>(graph.verticesNum);
		isEdgeAdd[i]= std::vector<bool>(graph.verticesNum);
	}

	for (int i = 0; i < graph.edgeNum; i++) { //заполнение матрицы смежности
		edgeMatrix[graph.edges[i].dot1][graph.edges[i].dot2] = true;
		edgeMatrix[graph.edges[i].dot2][graph.edges[i].dot1] = true;
	}

	for (int i = 0; i < graph.verticesNum; i++) { //заполнение секторов состоящих из двух вершин
		sectorMaxEdgeNum[i][i] = 0;
		sectorMaxEdgeNum[i][(i + 1) % graph.verticesNum] = edgeMatrix[enumerateVertices[i]][enumerateVertices[(i + 1) % graph.verticesNum]];
		isSectorWayExist[i][(i + 1) % graph.verticesNum]= edgeMatrix[enumerateVertices[i]][enumerateVertices[(i + 1) % graph.verticesNum]];
		isEdgeAdd[i][(i + 1) % graph.verticesNum]= edgeMatrix[enumerateVertices[i]][enumerateVertices[(i + 1) % graph.verticesNum]];
	}

	for (int i = 2; i < graph.verticesNum; i++) {
		for (int j = 0; j < graph.verticesNum; j++) {
			bool isWayExist = false;
			int maxEdgeNum = 0;
			int edge = -1;
			bool isEdgeAdded = false;
			for (int k = 1; k < i; k++) {
				bool flag = false;
				bool isEdge = false;
				int curEdgeNum = sectorMaxEdgeNum[j][(j + k) % graph.verticesNum] + sectorMaxEdgeNum[(j + k) % graph.verticesNum][(j + i) % graph.verticesNum]; 
				//количество рёбер при пересечении графа
				if (edgeMatrix[enumerateVertices[j]][enumerateVertices[(j + i) % graph.verticesNum]] && //проверка на возможности добавления ребра
					!(isSectorWayExist[j][(j + k) % graph.verticesNum] && isSectorWayExist[(j + k) % graph.verticesNum][(j + i) % graph.verticesNum])) {
					curEdgeNum += 1;
					flag = true;
					isEdge = true;
				}
				if (isSectorWayExist[j][(j + k) % graph.verticesNum] && isSectorWayExist[(j + k) % graph.verticesNum][(j + i) % graph.verticesNum]) { // проверка на существование пути в получившемся графе
					flag = true;
				}
				if (curEdgeNum > maxEdgeNum) { // выбор максимального из получвшихся подграфов в секторе
					maxEdgeNum = curEdgeNum;
					edge = j * graph.verticesNum + ((j + k) % graph.verticesNum);
					isWayExist = flag;
					isEdgeAdded = isEdge;
				}
			}
			sectorMaxEdgeNum[j][(j + i) % graph.verticesNum] = maxEdgeNum;
			sectorPointUnion[j][(j + i) % graph.verticesNum] = edge;
			isSectorWayExist[j][(j + i) % graph.verticesNum] = isWayExist;
			isEdgeAdd[j][(j + i) % graph.verticesNum] = isEdgeAdded;
		}
	}
	// обратный ход для восстановления последовательности рёбер
	int par1=0;
	int par2=0;
	int max = 0;
	for (int i = 0; i < graph.verticesNum; i++) {
		for (int j = 0; j < graph.verticesNum; j++) {
			if (sectorMaxEdgeNum[i][j] > max) {
				par1 = i;
				par2 = j;
				max = sectorMaxEdgeNum[i][j];
			}
		}
	}
	std::vector<Edge> edges= std::vector<Edge>();
	std::vector<int> queue = std::vector<int>();
	std::vector<int> queue2 = std::vector<int>();
	queue.push_back(par1);
	queue2.push_back(par2);
	int it = 0;
	while (it < queue.size()) {
		par1 = queue[it];
		par2 = queue2[it];
		if (isEdgeAdd[par1][par2]) {
			edges.push_back(Edge(enumerateVertices[par1], enumerateVertices[par2]));
		}
		if (sectorPointUnion[par1][par2] != -1) {
				int edge1 = sectorPointUnion[par1][par2] / graph.verticesNum;
				int edge2 = sectorPointUnion[par1][par2] % graph.verticesNum;
				queue.push_back(par1);
				queue2.push_back(edge2);
				queue.push_back(edge2);
				queue2.push_back(par2);
				}
		it++;
		}
	return Graph(graph.verticesNum, edges.size(), graph.vertices, edges);
}
bool isInOneCycle(int dot1, int dot2, int dot3, int dot4,std::vector<std::vector<int>>& listOfEdges) {
	std::vector<bool> isVisited(listOfEdges.size()); // Метка вершины
	std::vector<int> queue; // очередь для обхода в ширину
	std::vector<int> prev(listOfEdges.size(),-1); // предыдущая вершина
	queue.push_back(dot1);
	isVisited[dot1] = true;
	int j = 0;
	while (queue.size() > j) {
		int vert = queue[j];
		if (vert == dot2) {
			break;
		}
		for (int i = 0; i < listOfEdges[vert].size();i++) {
			if (!isVisited[listOfEdges[vert][i]]) {
				queue.push_back(listOfEdges[vert][i]);
				isVisited[listOfEdges[vert][i]] = true;
				prev[listOfEdges[vert][i]] = vert;
			}
		}
		j++;
	}
	int prevVert = prev[dot2];
	int curVert = dot2;
	while (prevVert != -1) { // проверка содержит ли путь от dot1 до dot2 ребро (dot3,dot4)
		if ((prevVert == dot3 && curVert == dot4) || (curVert == dot3 && prevVert == dot4)) {
			return true;
		}
		curVert = prevVert;
		prevVert = prev[prevVert];
	}
	return false;
}
Graph matroidAlgorithm(Graph graph) {
	DSU intersect(graph.edgeNum); // СНП для проверки, пересекают ли рёбра друг друга
	for (int i = 0; i < graph.edgeNum; i++) {
		if (intersect.getParrent(i) == i) {
			for (int j = 0; j < graph.crossings[i].size(); j++) {
				intersect.connect(i, graph.crossings[i][j]);
			}
		}
	}
	DSU graphCycle(graph.verticesNum); //СНП для недопущения циклов при построении начального множества
	std::vector<bool> intersectEdges(graph.edgeNum); // список для недопущения пересечений при построении начального множества
	std::vector<bool> resultEdges(graph.edgeNum); //множество J
	for (int i = 0; i < graph.edgeNum; i++) {
		if (!intersectEdges[intersect.getParrent(i)]) {
			if (graphCycle.connect(graph.edges[i].dot1, graph.edges[i].dot2)) {
				resultEdges[i] = true;
				intersectEdges[intersect.getParrent(i)] = true;
			}
		}
	}
	bool flag = true; // пометка для выхода из цикла
	while (flag) {
		flag = false;
		DSU graphCycle(graph.verticesNum); // построение СНП для проверки на циклы
		std::vector<std::vector<int>> listOfEdges(graph.verticesNum); // список переходов по рёбрам
		std::vector<bool> intersectEdges(graph.edgeNum); //список для проверки пересечения рёбер
		for (int i = 0; i < graph.edgeNum; i++) { // заполнение СНП, списка переходов, списка пересечений
			if (resultEdges[i]) {
				graphCycle.connect(graph.edges[i].dot1, graph.edges[i].dot2);
				intersectEdges[intersect.getParrent(i)] = true;
				listOfEdges[graph.edges[i].dot1].push_back(graph.edges[i].dot2);
				listOfEdges[graph.edges[i].dot2].push_back(graph.edges[i].dot1);
			}
		}
		std::vector<std::vector<int>> matroidEdges(graph.edgeNum); //рёбра в графе замен
		std::vector<bool> cycleMatroid(graph.edgeNum); //рёбра, которые можно добавить без появления циклов, множество X1
		std::vector<bool> intersectMatroid(graph.edgeNum); // рёбра, которые можно добавить без появления пересечений, множество X2
		for (int i = 0; i < graph.edgeNum; i++) { //построение ребёр графа пересечений и заполнение множеств X1,X2
			if (resultEdges[i]) {
				for (int j = 0; j < graph.edgeNum; j++) {
					if (!resultEdges[j]) { //добавление рёбер из J в S/J
						if (graphCycle.checkEdge(graph.edges[j].dot1, graph.edges[j].dot2)) {
							if (isInOneCycle(graph.edges[j].dot1, graph.edges[j].dot2, graph.edges[i].dot1, graph.edges[i].dot2,listOfEdges)) { //проверка на принадлежность одному циклу
								matroidEdges[i].push_back(j);
							}
						}
						else {
							matroidEdges[i].push_back(j);
						}
					}
				}
			}
			else {
				if (!graphCycle.checkEdge(graph.edges[i].dot1, graph.edges[i].dot2)) { //добавление ребра в X1
					cycleMatroid[i] = true;
				}
				if (!intersectEdges[intersect.getParrent(i)]) { //добалвение ребра в X2
					intersectMatroid[i] = true;
				}
				for (int j = 0; j < graph.edgeNum; j++) {// добавление рёбер из S/J в J
					if (resultEdges[j]) {
						if (intersect.checkEdge(i, j)) {
							matroidEdges[i].push_back(j); //замена ребра из пересечения на другое ребро из этого же пересечения
						}
						else {
							if (!intersectEdges[intersect.getParrent(j)]) {
								matroidEdges[i].push_back(j); //замена ребра из пересечения на ребро из пересечения, не содержащего ребра в J 
							}
						}
					}
				}
			}
		}
		std::vector<bool> isVisited(graph.edgeNum); // Метка рёбер
		std::vector<int> queue; // очередь для обработки рёбер обходом в ширину
		std::vector<int> prev(graph.edgeNum,-1); // предыдущее ребро
		for (int i = 0; i < graph.edgeNum; i++) {
			if (cycleMatroid[i]) {
				queue.push_back(i); //добавление ребёр из X1 в очередь и их пометка
				isVisited[i] = true;
			}
		}
		int j = 0;
		while (queue.size() > j) {
			int edge = queue[j];
			if (intersectMatroid[edge]) {
				flag = true; //пометка о том, что необоходимо продолжить итерации
				while (edge != -1) {
					resultEdges[edge] = (resultEdges[edge] != true); //изменение множества J
					edge = prev[edge];
				}
				break;
			}
			for (int i = 0; i < matroidEdges[edge].size(); i++) { // Добавление рёбер в очередь
				if (!isVisited[matroidEdges[edge][i]]) {
					queue.push_back(matroidEdges[edge][i]);
					prev[matroidEdges[edge][i]] = edge;
					isVisited[matroidEdges[edge][i]] = true;
				}
			}
			j++;
		}
	}
	std::vector<int> result; 
	for (int i = 0; i < graph.edgeNum; i++) { //запись из списка включений во множество рёбер
		if (resultEdges[i]) {
			result.push_back(i);
		}
	}
	return graph.getSubgraph(result); //возвращение подграфа
}
int main() {
	srand(time(NULL));
	int n = 1;
	double avg = 0;
for (int i = 0; i < n; i++) {
	Graph graph = createRandomCircleGraph(5, 5, 2);
		auto begin = std::chrono::steady_clock::now();
		graph.getCrossing();
		int crossNum = 0;
		for (int i = 0; i < graph.edgeNum; i++) {
			crossNum += graph.crossings[i].size();
		}
		std::cout << "Crossings: " << crossNum << '\n';
		//Graph subGraph = baseAlgorithm(graph);
		Graph subGraph2 = matroidAlgorithm(graph);
		std::cout << "Iteration: " << std::to_string(i + 1) << '\n';
		auto end = std::chrono::steady_clock::now();
		graph.saveToTextFile("graph" + std::to_string(i + 1) + ".txt");
		//subGraph.saveToTextFile("graph1" + std::to_string(i + 1) + ".txt");
		//subGraph2.saveToTextFile("graph2" + std::to_string(i + 1) + ".txt");
		drawGraph("graph" + std::to_string(i + 1) + ".txt", "graph" + std::to_string(i + 1) + ".png", "black");
		//drawGraph("graph1" + std::to_string(i + 1) + ".txt", "graph1" + std::to_string(i + 1) + ".png", "black");
		//drawGraph("graph2" + std::to_string(i + 1) + ".txt", "graph2" + std::to_string(i + 1) + ".png", "black");
		avg += (std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() * 1.0) / n;
	}
	std::cout<<"Average time:" << avg << "ms\n";
}