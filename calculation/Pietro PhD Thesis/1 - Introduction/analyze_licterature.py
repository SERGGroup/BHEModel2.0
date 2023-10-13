# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_code.constants import CALCULATION_FOLDER, os
from collections.abc import Iterable
from abc import ABC, abstractmethod
from pyvis.network import Network
from collections import deque
import pandas as pd
import numpy as np


class BaseTree(ABC):

    def __init__(self):

        self.left = None
        self.right = None
        self.__elements = None
        self.__tmp_counter = 0.

    @classmethod
    def init_new_element(cls, new_value):

        new_class = cls()
        new_class.initialize_class_values(new_value)
        return new_class

    def append_new_value(self, new_value):

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

    def convert_in_list(self):

        return_list = list()
        @self.while_visiting()
        def counting_function(current_node):
            return_list.append(current_node)

        counting_function()
        return return_list

    def count(self, bool_func):

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

    def append_author_list(self, author_list):

        if self.surname is None:

            self.initialize_class_values(author_list[0])
            author_list = author_list[1:]

        for author in author_list:
            self.append_new_value(author)

    def __init__(self):
        super().__init__()

        self.names = list()
        self.surname = None

    def initialize_class_values(self, new_value):

        if len(new_value) == 1:

            self.names = [""]
            self.surname = new_value[0]

        else:

            self.names = [new_value[0]]
            self.surname = new_value[1]

    def update_if_append_equal(self, new_value):
        if new_value[0] not in self.names:
            self.names.append(new_value[0])

    def convert_in_list(self):

        tmp_list = super().convert_in_list()

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

            return self.surname == other

        return self.surname == other.surname

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

        return "{}".format(self.surname)


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

    def generate_author_network(self, tmp_author_tree=None):

        if tmp_author_tree is None:

            tmp_author_tree = self.generate_author_tree()

        tmp_author_list = tmp_author_tree.convert_in_list()
        new_network = Network()
        new_network.add_nodes(tmp_author_list)

        @self.while_visiting()
        def append_author_list(current_node: ArticleTree):

            tmp_authors = current_node.authors
            tmp_n_authors = len(tmp_authors)

            if tmp_n_authors > 1:

                for i in range(tmp_n_authors):

                    source_author = tmp_author_tree.search(tmp_authors[i])

                    for j in range(tmp_n_authors - i - 1):

                        to_author = tmp_author_tree.search(tmp_authors[j + i + 1])
                        new_network.add_edge(str(source_author), str(to_author))

        append_author_list()
        return new_network

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

        if self.title == other.title:
            return self.citations > self.citations

        return self.title > other.title

    def __lt__(self, other):

        # enables comparison
        # self < other
        if self.title == other.title:
            return self.citations < self.citations

        return self.title < other.title

    def __le__(self, other):

        return not self.__gt__(other)

    def __ge__(self, other):

        return not self.__lt__(other)

    def __str__(self):

        return "{}".format(self.title)


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

# %%-------------------------------------   PERFORM CALCULATIONS                -------------------------------------> #
net = article_tree.generate_author_network()
html = net.show(os.path.join(resources_folder, "Authors Networks.html"), notebook=False)