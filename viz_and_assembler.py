from Bio import SeqIO
from graphviz import Digraph
import argparse
import pydot

class Vertex:

    def __init__(self, seq):
        self.seq = seq
        self.vertex_coverage = 1
        self.in_edges = {}
        self.out_edges = {}
    def increase_vertex_coverage(self):
        self.vertex_coverage += 1

class Edge:

    def __init__(self, current_kmer, next_kmer):
        self.seq = current_kmer + next_kmer[-1]
        self.in_vertices = {}
        self.out_vertices = {}
        self.edge_coverage = 0
    def increase_edge_coverage(self):
        self.edge_coverage += 1
    def calculation_edge_coverage(self, prev_vertex_cov, next_vertext_cov):
        self.edge_coverage = (prev_vertex_cov + next_vertext_cov)/2

class Graph:

    def __init__(self, k):
        self.vertices = {}
        self.k = k
    def add_read(self, read):

        if len(read) < self.k:
            return

        # first k-mer
        kmer = read[:self.k]
        if kmer in self.vertices:
            self.vertices[kmer].increase_vertex_coverage()
        else:
            self.vertices[kmer] = Vertex(kmer)

        # next k-mer iterations:
        for i in range(1, len(read)-self.k+1, 1):
            next_kmer = read[i:i+self.k]
            if next_kmer in self.vertices:
                self.vertices[next_kmer].increase_vertex_coverage()
            else:
                self.vertices[next_kmer] = Vertex(next_kmer)

            # add new edge
            new_edge = Edge(kmer, next_kmer)
            # add vertices
            self.vertices[next_kmer].in_edges[kmer] = [new_edge]
            self.vertices[kmer].out_edges[next_kmer] = [new_edge]
            kmer = next_kmer
    def coverage_calculating(self):
        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].out_edges.keys():
                self.vertices[current_vertex].out_edges[next_vertex][0]\
                      .calculation_edge_coverage(self.vertices[current_vertex].vertex_coverage,
                      self.vertices[next_vertex].vertex_coverage)
    def visualize_graph(self, arg):
        dot = Digraph(comment='Assemble')
        if arg == 'full':
            for k, v in self.vertices.items():
                dot.node(k, label='{}'.format(k))
                for kk, vv in v.out_edges.items():
                    dot.edge(k, kk, label='{}'.format(vv[0].seq))
        else:
            for k, v in self.vertices.items():
                dot.node(k, label='cov={}'.format(v.vertex_coverage))
                for kk, vv in v.out_edges.items():
                    dot.edge(k, kk, label='coverage:{} length:{}'.format(vv[0].edge_coverage,len(vv[0].seq)))
        dot.format = 'svg'
        dot.render('de_bjuijn.svg', view=True)

    def launch_assembler(self):
        begin = []
        for i in list(self.vertices.keys()):
            if (len(self.vertices[i].out_edges.keys()) == 1) & (len(self.vertices[i].in_edges.keys()) == 1):
                begin.append(i)
        for i in begin:
            prev = list(self.vertices[i].in_edges.keys())[0]
            next = list(self.vertices[i].out_edges.keys())[0]
            self.vertices[prev].out_edges[next] = [Edge(prev, next)]
            self.vertices[next].in_edges[prev] = [Edge(prev, next)]
            self.vertices[next].in_edges[prev][0].seq = self.vertices[i].in_edges[prev][0].seq + self.vertices[next].in_edges[i][0].seq[self.k:]
            self.vertices[prev].out_edges[next][0].seq = self.vertices[i].in_edges[prev][0].seq + self.vertices[next].in_edges[i][0].seq[self.k:]
            cov1 = self.vertices[i].out_edges[next][0].edge_coverage
            cov2 = self.vertices[i].in_edges[prev][0].edge_coverage
            cov_vertex = self.vertices[i].vertex_coverage
            len1 = len(self.vertices[i].out_edges[next][0].seq)
            len2 = len(self.vertices[i].in_edges[prev][0].seq)
            len3 = len(self.vertices[prev].out_edges[next][0].seq)
            edge_coverage_2 = (cov1 * (len1 - self.k + 1) + cov2 * (len2 - self.k + 1) - cov_vertex) / (len3 - self.k + 1)
            self.vertices[prev].out_edges[next][0].coverage = edge_coverage_2
            self.vertices[next].in_edges[prev][0].coverage = edge_coverage_2
            del self.vertices[prev].out_edges[i]
            del self.vertices[next].in_edges[i]
            del self.vertices[i]
        print('Assemble finished')

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',type=str)
    parser.add_argument('-k', '--kmer', help='k-mer length', default=15, type=int)
    parser.add_argument('-t', '--type', help='visualization type: full/short', default='full', type=str)
    parser.add_argument('-s', '--strand', help='fw|bw', default='fw', type=str)

    args = parser.parse_args()
    my_graph = Graph(args.kmer)
    if args.strand == 'fw':
        with open(args.input) as f:
            for record in SeqIO.parse(f, 'fasta'):
                my_graph.add_read(str(record.seq))
    else:
        with open(args.input) as f:
            for record in SeqIO.parse(f, 'fasta'):
                my_graph.add_read(str(record.reverse_complement().seq))

    #my_graph.coverage_calculating()
    my_graph.visualize_graph(arg=args.type)
    #my_graph.launch_assembler()