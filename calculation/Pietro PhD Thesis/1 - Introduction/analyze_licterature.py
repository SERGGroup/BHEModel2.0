# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.constants import CALCULATION_FOLDER, os
from collections.abc import Iterable
from abc import ABC, abstractmethod
from pyvis.network import Network
from collections import deque
import pandas as pd
import numpy as np


# %%-------------------------------------   DEFINE CLASSES                      -------------------------------------> #
class BaseTree(ABC):

    def __init__(self):

        self.left = None
        self.right = None
        self.init_value = None
        self.__elements = None
        self.__tmp_counter = 0.

    @classmethod
    def init_new_element(cls, new_value):

        new_class = cls()
        new_class.init_value = new_value
        new_class.initialize_class_values(new_value)
        return new_class

    def append_new_value(self, new_value):

        if self.init_value is None:

            self.init_value = new_value
            self.initialize_class_values(new_value)

        else:

            new_point = type(self).init_new_element(new_value)
            curr_point = self

            while curr_point is not None:

                curr_point.__elements = None

                if curr_point == new_point:

                    curr_point.update_if_append_equal(new_value)

                    break

                elif curr_point > new_point:

                    if curr_point.left is None:

                        curr_point.left = new_point
                        break

                    else:

                        curr_point = curr_point.left

                else:

                    if curr_point.right is None:

                        curr_point.right = new_point
                        break

                    else:
                        curr_point = curr_point.right

    def search(self, value):

        curr_point = self

        while curr_point is not None:

            if curr_point == value:

                return curr_point

            elif curr_point > value:

                curr_point = curr_point.left

            else:

                curr_point = curr_point.right

    def while_visiting(self, run_position="center"):

        def wrapper_function(func):

            def visiting_function(*args, **kwargs):

                stack = deque()
                current_node = self

                while stack or current_node:

                    if current_node:

                        if run_position not in ["left", "right"]:
                            func(current_node, *args, **kwargs)

                        stack.append(current_node)
                        current_node = current_node.left

                    else:

                        current_node = stack.pop()
                        if run_position == "left":
                            func(*args, **kwargs)

                        current_node = current_node.right
                        if run_position == "right":
                            func(*args, **kwargs)

            return visiting_function

        return wrapper_function

    def list_values(self, bool_func=None):

        if bool_func is None:
            def bool_func(val):
                return True

        return_list = list()
        @self.while_visiting()
        def listing_function(current_node):
            if bool_func(current_node):
                return_list.append(current_node)

        listing_function()
        return return_list

    def count(self, bool_func=None):

        if bool_func is None:
            def bool_func(val):
                return True

        self.__tmp_counter = 0
        @self.while_visiting()
        def counting_function(current_node):
            if bool_func(current_node):
                self.__tmp_counter += 1

        counting_function()
        return self.__tmp_counter

    @property
    def elements(self):

        if self.__elements is None:

            self.__elements = 1

            if self.left is not None:

                self.__elements += self.left.elements

            if self.right is not None:

                self.__elements += self.right.elements

        return self.__elements

    @abstractmethod
    def initialize_class_values(self, new_value):
        pass

    @abstractmethod
    def update_if_append_equal(self, new_value):
        pass


class AuthorTree(BaseTree):

    def append_author_list(self, author_list_in):

        if self.surname is None:

            self.initialize_class_values(author_list_in[0])
            author_list_in = author_list_in[1:]

        for author in author_list_in:
            self.append_new_value(author)

    def __init__(self):
        super().__init__()

        self.names = list()
        self.surname = None

    @property
    def name(self):
        if len(self.names) > 0:
            return self.names[0]
        else:
            return ""

    def initialize_class_values(self, new_value):

        self.init_value = new_value

        if len(new_value) == 1:

            self.names = [""]
            self.surname = new_value[0].strip("...")

        else:

            self.names = [new_value[0].strip("...")]
            self.surname = new_value[1].strip("...")

    def update_if_append_equal(self, new_value):
        if new_value[0] not in self.names:
            self.names.append(new_value[0])

    def list_values_str(self, bool_func=None):

        tmp_list = self.list_values(bool_func)

        new_list = list()
        for elem in tmp_list:

            new_list.append(elem.surname)

        return new_list

    def search(self, value):

        if isinstance(value, Iterable):

            tmp_value = AuthorTree()
            tmp_value.initialize_class_values(value)
            value = tmp_value

        return super().search(value)

    def __eq__(self, other):

        if type(other) == str:

            return str(self) == other

        return str(self) == str(other)

    def __gt__(self, other):
        # enables comparison
        # self > other

        if type(other) == str:
            return self.surname > other

        return self.surname > other.surname

    def __lt__(self, other):
        # enables comparison
        # self < other
        if type(other) == str:

            return self.surname < other

        return self.surname < other.surname

    def __le__(self, other):
        return not self.__gt__(other)

    def __ge__(self, other):
        return not self.__lt__(other)

    def __str__(self):

        return_str = "{}".format(self.init_value[0])

        for i in range(len(self.init_value) - 1):

            return_str += ", {}".format(self.init_value[i + 1])

        return return_str


class ArticleTree(BaseTree):

    @classmethod
    def init_from_df(cls, df_in):

        new_element = cls()
        new_element.initialize_class_values(df_in.iloc[0])
        new_element.append_new_values(df_in.iloc[1:])

        return new_element

    def append_new_values(self, new_values):

        for i in range(len(new_values.Title)):

            self.append_new_value(new_values.iloc[i])

    def __init__(self):

        super().__init__()

        self.year = None
        self.title = None
        self.citations = None
        self.citation_author = None

        self.n_authors = 0.
        self.__authors = list()
        self.other_authors = False

        self.cite_brown = False
        self.cite_pruess = False

    def initialize_class_values(self, new_value):

        self.year = new_value.Year
        self.title = new_value.Title
        self.citations = new_value.Cites
        self.citation_author = new_value.CitesPerAuthor

        self.n_authors = new_value.AuthorCount
        self.authors = new_value.Authors

        if new_value.IsBrown:
            self.cite_brown = True
        else:
            self.cite_pruess = True

    def update_if_append_equal(self, new_value):

        if new_value.IsBrown:
            self.cite_brown = True

        else:
            self.cite_pruess = True

    def generate_author_tree(self):

        new_author_tree = AuthorTree()

        @self.while_visiting()
        def append_author_list(current_node: ArticleTree):

            new_author_tree.append_author_list(current_node.authors)

        append_author_list()
        return new_author_tree

    @property
    def authors(self):

        return self.__authors

    @authors.setter
    def authors(self, autor_raw):

        self.__authors = list()
        self.__init_author_list(autor_raw)

    def __init_author_list(self, autor_raw):

        for author in autor_raw.split(", "):

            if author == '...':

                self.other_authors = True

            else:

                self.authors.append(author.split())

    def __eq__(self, other):

        return self.title == other.title and self.citations == other.citations

    def __gt__(self, other):

        # enables comparison
        # self > other

        if self.citations == other.citations:
            return self.title > self.title

        return self.citations < other.citations

    def __lt__(self, other):

        # enables comparison
        # self < other
        if self.citations == other.citations:
            return self.title < self.title

        return self.citations > other.citations

    def __le__(self, other):

        return not self.__gt__(other)

    def __ge__(self, other):

        return not self.__lt__(other)

    def __str__(self):

        return "{}".format(self.title)


class AuthorsNetwork:

    class AuthorNode(AuthorTree):

        def __init__(self):

            super().__init__()
            self.connections = list()

        @classmethod
        def __init_from_author_tree(cls, input_tree: AuthorTree):

            new_cls = cls()
            new_cls.append_new_value(input_tree.init_value)

        @property
        def citation_score(self):

            score = 0.

            for connection in self.connections:

                score += connection.n_cit_cum

            return score

    class AuthorConnection(BaseTree):

        def initialize_class_values(self, new_value):

            if new_value[0] > new_value[1]:

                self.from_author = new_value[0]
                self.to_author = new_value[1]

            else:

                self.from_author = new_value[1]
                self.to_author = new_value[0]

            new_value[0].connections.append(self)
            new_value[1].connections.append(self)

            self.n_cit_cum = new_value[2]
            self.n_article_cum = 1

        def update_if_append_equal(self, new_value):

            self.n_article_cum += 1
            self.n_cit_cum += new_value[2]

        def __init__(self):

            super().__init__()
            self.from_author = None
            self.to_author = None
            self.n_article_cum = 0.
            self.n_cit_cum = 0.

        def __eq__(self, other):

            return (self.from_author == other.from_author) and (self.to_author == other.to_author)

        def __gt__(self, other):

            # enables comparison
            # self > other
            if self.from_author == other.from_author:
                return self.to_author > other.to_author
            return self.from_author > other.from_author

        def __lt__(self, other):

            # enables comparison
            # self < other
            if self.from_author == other.from_author:
                return self.to_author < other.to_author
            return self.from_author < other.from_author

        def __le__(self, other):
            return not self.__gt__(other)

        def __ge__(self, other):
            return not self.__lt__(other)

    def __init__(self, author_list_in, article_list_in):

        self.author_list = author_list_in
        self.article_list = article_list_in

        self.author_tree = self.AuthorNode()
        self.connection_tree = self.AuthorConnection()
        self.plt_network = None

        self.__init_nodes()
        self.__init_connections()

    def __init_nodes(self):

        for author in self.author_list:

            self.author_tree.append_new_value(author.init_value)

    def __init_connections(self):

        for article in self.article_list:

            tmp_authors = article.authors
            tmp_n_authors = len(tmp_authors)

            if tmp_n_authors > 1:

                for i in range(tmp_n_authors):

                    source_author = self.author_tree.search(tmp_authors[i])

                    if source_author is not None:

                        for j in range(tmp_n_authors - i - 1):

                            to_author = self.author_tree.search(tmp_authors[j + i + 1])

                            if to_author is not None:

                                self.connection_tree.append_new_value([

                                    source_author,
                                    to_author,
                                    article.citation_author

                                ])

    def init_network(self, min_connections=0):

        self.plt_network = Network()

        @self.author_tree.while_visiting()
        def append_author_list(current_node):

            if len(current_node.connections) > min_connections:

                self.plt_network.add_node(

                    str(current_node),
                    size=np.log(current_node.citation_score)*10,
                    label=current_node.surname

                )

        @self.connection_tree.while_visiting()
        def append_connections(current_node):

            if (len(current_node.from_author.connections) > min_connections) and \
                    (len(current_node.to_author.connections) > min_connections):

                self.plt_network.add_edge(

                    str(current_node.from_author),
                    str(current_node.to_author),
                    width=np.log(current_node.n_cit_cum)*5

                )

        append_author_list()
        append_connections()

    def show(self):

        if self.plt_network is None:
            self.init_network()

        self.plt_network.show(

            os.path.join(resources_folder, "Authors Networks.html"),
            notebook=False

        )


# %%-------------------------------------   IMPORT EXCEL FILE                   -------------------------------------> #
resources_folder = os.path.join(

    CALCULATION_FOLDER, "Pietro PhD Thesis",
    "1 - Introduction", "res"

)

excel_path = os.path.join(resources_folder, "Licterature Review - Data.xlsx")

brown_data = pd.read_excel(excel_path, sheet_name="DATA")
pruess_data = pd.read_excel(excel_path, sheet_name="DATA Pruess")

article_tree = ArticleTree.init_from_df(brown_data)
article_tree.append_new_values(pruess_data)
author_tree = article_tree.generate_author_tree()


# %%-------------------------------------   PERFORM CALCULATIONS                -------------------------------------> #
@article_tree.count
def with_year(current_node):
    return not np.isnan(current_node.year)


@article_tree.count
def cite_both(current_node):
    return current_node.cite_brown and current_node.cite_pruess


@article_tree.count
def cite_brown(current_node):
    return current_node.cite_brown


@article_tree.count
def cite_pruess(current_node):
    return current_node.cite_pruess


print(with_year)
print(cite_both)
print(cite_brown - cite_both)
print(cite_pruess - cite_both)


# %%-------------------------------------   INIT NETWORK                        -------------------------------------> #
@article_tree.list_values
def article_list(current_node):
    return True


@author_tree.list_values
def author_list(current_node):
    return len(current_node.surname) > 3


network = AuthorsNetwork(

    author_list,
    article_list

)

network.init_network(min_connections=2)
network.show()


# %%-------------------------------------   PRINT ARTICLES                      -------------------------------------> #
min_citations = -1
@author_tree.list_values
def author_saar_list(current_node):
    return current_node.surname in ["Saar", "Adams", "Randolph", "Bielicki"]

@article_tree.list_values
def article_saar_list(current_node):

    for author_saar in author_saar_list:

        if author_saar.init_value in current_node.authors:

            if current_node.citations > min_citations:
                return True

    return False

# %%-------------------------------------   PRINT ARTICLES                      -------------------------------------> #
min_citations = -1
@author_tree.list_values
def author_yao_list(current_node):
    return current_node.init_value in [["F", "Sun"]]

@article_tree.list_values
def article_yao_list(current_node):

    for author_yao in author_yao_list:

        if author_yao.init_value in current_node.authors:

            if current_node.citations > min_citations:
                return True

    return False

citations = 0.
for article in article_yao_list:

    citations += article.citations

print(citations)